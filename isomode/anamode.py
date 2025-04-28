import tempfile
import os
import numpy as np
from abipy.abilab import abiopen
from ase.io import read, write
from isomode.pydistort import isocif, isodistort, view_distort, logger
from isomode.phonmode import phonon_distort_generator, displacement_cart_to_evec
from isomode.isodistort_parser import Isomode
from isomode.frozen_mode import get_supercell, distorted_cell

class LabelPhbst(object):
    """
    A class for labeling phonon modes from ABINIT PHBST calculations using ISODISTORT.

    This class analyzes phonon modes from ABINIT calculations, matches them with
    symmetry-adapted modes from ISODISTORT, and provides labeling information.

    Attributes:
        fname (str): Path to the PHBST.nc file
        tmpdir (str): Temporary directory for intermediate files
        sc_mat (ndarray): Supercell matrix used for analysis
        qdict (dict): Dictionary of q-points to analyze
        parent_atoms (ase.Atoms): Parent structure
        mode_info (dict): Collected mode information (frequencies, labels, symmetries)
    """
    def __init__(self,
                 fname,
                 tmpdir=None,
                 sc_mat=np.eye(3) * 2,
                 qdict={
                     'Gamma': [0.0, 0.0, 0.0],
                 }):
        self.fname = fname
        self.sc_mat = sc_mat
        self.qdict = qdict
        self._read_file()

        if tmpdir is None:
            self._tmpdir = tempfile.mkdtemp()
        else:
            self._tmpdir = tmpdir
            if not os.path.exists(self._tmpdir):
                os.makedirs(self._tmpdir)

    def _read_file(self):
        """
        Read and parse the PHBST.nc file.

        Extracts:
        - Atomic structure
        - Q-points
        - Phonon displacements
        - Frequencies (converted to cm-1)
        - Atomic masses
        """
        abifile = abiopen(self.fname)
        self.parent_atoms = abifile.structure.to_ase_atoms()
        self.parent_atoms.set_pbc([True, True, True])

        self.qpoints = abifile.phbands.qpoints.frac_coords
        self.displacements = abifile.phbands.phdispl_cart
        self.freqs=abifile.phbands.phfreqs*8056.6
        self.masses = self.parent_atoms.get_masses()

        self.idq = {}
        self.pdisps = {}
        self.edisps = {}
        self.efreqs = {}
        for name, qpt in self.qdict.items():
            for i, q in enumerate(self.qpoints):
                if np.allclose(qpt, q):
                    self.idq[name] = i
                    self.edisps[name] = self.displacements[i]
                    self.efreqs[name] = self.freqs[i]
            if name not in self.idq:
                print("WARNING: The qpoint %s: %s not found in phonon" % (name, qpt))
        self.dcell = distorted_cell(
            self.parent_atoms, supercell_matrix=self.sc_mat)

    def _write_cif(self):
        """
        Write structure to CIF and generate symmetrized CIF file.

        Creates two files in tmpdir:
        - parent.cif: Original structure
        - parent_sym.cif: Symmetrized structure

        The symmetrized CIF is used for ISODISTORT analysis.
        """
        parent_cif = os.path.join(self._tmpdir, 'parent.cif')
        write(parent_cif, self.parent_atoms)
        # convert highsym_fname
        isosym = isocif(parent_cif)
        isosym.upload_cif()
        isosym.findsym()
        self.parent_sym_cif = os.path.join(self._tmpdir, 'parent_sym.cif')
        isosym.save_cif(fname=self.parent_sym_cif)

    def _write_supercell_cif(self):
        """
        Generate and write supercell structure to CIF file.

        Creates:
        - supercell.cif: Structure expanded by sc_mat

        The supercell is needed for mode visualization and analysis.
        """
        self.sc_atoms = get_supercell(
            self.parent_atoms, self.sc_mat, symprec=1e-5)

        self.sc_atoms.set_pbc([True, True, True])
        self.supercell_cif = os.path.join(self._tmpdir, 'supercell.cif')
        write(self.supercell_cif, self.sc_atoms)

    def _get_mode_details(self):
        """
        Run ISODISTORT to get mode symmetry information.

        Creates:
        - mode_details.txt: File containing ISODISTORT mode definitions

        This connects the phonon modes to irreducible representations.
        """
        d = isodistort(
            parent_cif=self.parent_sym_cif, distorted_cif=self.supercell_cif)
        self.mode_fname = os.path.join(self._tmpdir, 'mode_details.txt')
        d.get_mode_details(save_fname=self.mode_fname)

    def get_modes(self):
        """
        Parse ISODISTORT mode definitions file.

        Populates:
        - mode_definitions: Dictionary of mode information including:
          * symmetry labels
          * displacement patterns
          * q-point information
        """
        m = Isomode(self.mode_fname)
        m.read_mode_definitions()
        #self.mode_displacements=m.get_mode_displacement()
        self.mode_definitions = m.mode_definitions

    def displ_cart_to_supercell_disp(self, displ_cart, qpt):
        """
        Convert phonon displacement to supercell displacement pattern.

        Args:
            displ_cart (ndarray): Cartesian displacements (3N vector)
            qpt (list): Q-point in fractional coordinates

        Returns:
            ndarray: Displacement pattern in supercell
        """
        evec = displacement_cart_to_evec(
            displ_cart,
            masses=self.masses/self.masses,
            scaled_positions=self.parent_atoms.get_scaled_positions(),
            qpoint=qpt,
            add_phase=True)
        disp = self.dcell._get_displacements(
            evec, qpt, amplitude=0.5, argument=0)
        return disp

    def projection(self):
        """
        Project phonon modes onto ISODISTORT symmetry modes.

        Returns:
            dict: Mapping of (qname, mode_index) to:
                  * frequency (cm-1)
                  * symmetry label
                  * irreducible representation
        """
        logger.info("Projecting phonon modes onto symmetry modes")
        ret = {}
        masses = self.sc_atoms.get_masses()
        for qname, d in self.edisps.items():
            qpt = self.qdict[qname]
            for imode, mode in enumerate(d):
                freq=self.efreqs[qname][imode]
                sc_disp = self.displ_cart_to_supercell_disp(mode,
                                                            qpt).flatten()
                #sc_evec = np.array(sc_disp) * np.sqrt(
                #    np.kron(masses, [1, 1, 1]))
                sc_evec = np.real(sc_disp) / np.linalg.norm(sc_disp)
                #print("================ imode: %s , freqs: %s ========== "%(imode,freq ) )
                #print(np.real(sc_disp).reshape((40,3)))
                #print(np.real(mode))
                for fullname, md in self.mode_definitions.items():
                    mdkpt = md['mode_info']['kpt']
                    mdlabel = md['mode_info']['label']
                    mdsym =md['mode_info']['symmetry']
                    mddisp = md['mode'].flatten()

                    #mdevec = np.array(mddisp) * np.sqrt(
                    #    np.kron(masses, [1, 1, 1]))
                    mdevec=mddisp
                    mdevec = np.real(mdevec) / np.linalg.norm(mdevec)

                    proj = np.abs(np.dot(sc_evec, mdevec))
                    #if np.allclose(mdkpt, qpt) and proj > 0.5:
                    if proj>0.5:
                        #print("md_evec", mdevec)
                        #print("sc_evec", sc_evec)
                        #print(fullname)
                        #print(mdevec.reshape((len(self.sc_atoms), 3)))
                        #print(mdlabel)
                        #print(proj)
                        ret[(qname, imode)] = {'freq': freq, 'label': mdlabel,  'symmetry': mdsym}
        self.mode_info=ret
        logger.info("Projection complete, Labeling phonon modes by projections")
        return ret

    def get_labels(self):
        """
        Perform complete labeling workflow.

        Steps:
        1. Write CIF files
        2. Run ISODISTORT analysis
        3. Get mode definitions
        4. Project phonons onto symmetry modes

        Returns:
            dict: Mode labeling information (same as projection())
        """
        self._write_cif()
        self._write_supercell_cif()
        self._get_mode_details()
        self.get_modes()
        return self.projection()

    def write_result(self, fname):
        """
        Write mode labeling results to file.

        Args:
            fname (str): Output filename

        File format:
        # qname   id  freq(cm-1)  label  symmetry
        Gamma     0   100.0       GM1+   A1g
        Gamma     1   150.0       GM2-   Eg
        ...
        """
        self.get_labels()
        with open(fname, 'w') as myfile:
            myfile.write("# %7s  %5s  %10s  %5s  %9s\n"%("qname",'id', 'freq(cm-1)', 'label','symmetry'))
            for key, val in self.mode_info.items():
                qname, imode =key
                freq=val['freq']
                label=val['label']
                sym=val['symmetry']
                myfile.write("%8s  %5d  %10.5e  %5s  %9s\n"%(qname, imode, freq, label, sym))

def label_phbst(fname_phbst, output_fname, tmpdir=None,
                         sc_mat=np.eye(3) * 2,
                         qdict={'Gamma': [0.0, 0.0, 0.0]}):
    """
    Combined function to initialize LabelPhbst and write results to a file in one step.
    
    This function simplifies the process of analyzing phonon modes by combining the
    initialization of LabelPhbst with writing the results to a file.
    
    Parameters:
    -----------
    fname_phbst : str
        Path to the ABINIT PHBST.nc file containing phonon data.
        Expected to be generated from a PHBST calculation.
        
    output_fname : str
        Path where the output results file will be written.
        File will contain mode frequencies, labels and symmetries.
        
    tmpdir : str, optional
        Temporary directory for intermediate files (default: None).
        If None, a temporary directory will be created automatically.
        
    sc_mat : numpy.ndarray, optional
        Supercell matrix used for mode analysis (default: 2x2x2).
        Should be a 3x3 numpy array specifying the supercell transformation.
        
    qdict : dict, optional
        Dictionary of q-points to analyze (default: {'Gamma': [0.0, 0.0, 0.0]}).
        Format: {'label': [q_x, q_y, q_z]} where label is a string identifier.
        
    Returns:
    --------
    None
    
    Examples:
    --------
    >>> # Basic usage with default parameters
    >>> label_phbst_and_write('phonons/run.abo_PHBST.nc', 'results.txt')
    
    >>> # Custom supercell and q-points
    >>> import numpy as np
    >>> sc = np.eye(3) * 3  # 3x3x3 supercell
    >>> qpts = {'Gamma': [0,0,0], 'X': [0.5,0,0]}
    >>> label_phbst_and_write('phonons.nc', 'modes.txt', sc_mat=sc, qdict=qpts)
    """
    l = LabelPhbst(fname_phbst, tmpdir=tmpdir, sc_mat=sc_mat, qdict=qdict)
    l.write_result(output_fname)

def main():
    """
    Command line interface for phonon mode labeling.
    
    Usage examples:
    python anamode.py input.nc output.txt
    python anamode.py input.nc output.txt --tmpdir ./temp
    python anamode.py input.nc output.txt --sc_mat 2 2 2 --qpoints Gamma 0 0 0 X 0.5 0 0
    """
    import argparse
    import numpy as np
    
    parser = argparse.ArgumentParser(description='Label phonon modes from ABINIT PHBST file')
    parser.add_argument('input', help='Input PHBST.nc file')
    parser.add_argument('output', help='Output results file')
    parser.add_argument('--tmpdir', help='Temporary directory (default: auto-create)')
    parser.add_argument('--sc_mat', nargs=3, type=int, default=[2,2,2],
                       help='Supercell matrix as 3 integers (default: 2 2 2)')
    parser.add_argument('--qpoints', nargs='+',
                       help='Q-points as alternating labels and coordinates (e.g. Gamma 0 0 0 X 0.5 0 0)')
    
    args = parser.parse_args()
    
    # Process qpoints if provided
    qdict = {'Gamma': [0.0, 0.0, 0.0]}  # default
    if args.qpoints:
        if len(args.qpoints) % 4 != 0:
            raise ValueError("Q-points must be in groups of 4 (label x y z)")
        qdict = {}
        for i in range(0, len(args.qpoints), 4):
            label = args.qpoints[i]
            coords = list(map(float, args.qpoints[i+1:i+4]))
            qdict[label] = coords
    
    # Create supercell matrix
    sc_mat = np.diag(args.sc_mat)
    
    label_phbst(
        args.input,
        args.output,
        tmpdir=args.tmpdir,
        sc_mat=sc_mat,
        qdict=qdict
    )

if __name__ == '__main__':
    main()
