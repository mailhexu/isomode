import tempfile
import os
import numpy as np
from abipy.abilab import abiopen
from ase.io import read, write
from isomode.pydistort import isocif, isodistort, view_distort
from isomode.phonmode import phonon_distort_generator, displacement_cart_to_evec
from isomode.isodistort_parser import Isomode
from isomode.frozen_mode import get_supercell, distorted_cell


class LabelDDB(object):
    def __init__(self,
                 fname,
                 tmpdir=None,
                 sc_mat=np.eye(3) * 2,
                 qdict={
                     'Gamma': [0.0, 0.0, 0.0],
                     'R': [0.5, 0.5, 0.5],
                     'M':[0.5,0.5,0],
                     'X':[0,0.5,0]
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
                print("The qpoint %s: %s not found in phonon" % (name, qpt))
        self.dcell = distorted_cell(
            self.parent_atoms, supercell_matrix=self.sc_mat)

    def _write_cif(self):
        """
        write structure to cif and the symmetrized cif file.
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
        self.sc_atoms = get_supercell(
            self.parent_atoms, self.sc_mat, symprec=1e-5)

        self.sc_atoms.set_pbc([True, True, True])
        self.supercell_cif = os.path.join(self._tmpdir, 'supercell.cif')
        write(self.supercell_cif, self.sc_atoms)

    def _get_mode_details(self):
        d = isodistort(
            parent_cif=self.parent_sym_cif, distorted_cif=self.supercell_cif)
        self.mode_fname = os.path.join(self._tmpdir, 'mode_details.txt')
        d.get_mode_details(save_fname=self.mode_fname)

    def get_modes(self):
        m = Isomode(self.mode_fname)
        m.read_mode_definitions()
        #self.mode_displacements=m.get_mode_displacement()
        self.mode_definitions = m.mode_definitions

    def displ_cart_to_supercell_disp(self, displ_cart, qpt):
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
                print("================ imode: %s , freqs: %s ========== "%(imode,freq ) )
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
                        print(fullname)
                        #print(mdevec.reshape((len(self.sc_atoms), 3)))
                        print(mdlabel)
                        print(proj)
                        ret[(qname, imode)] = {'freq': freq, 'label': mdlabel,  'symmetry': mdsym}
        self.mode_info=ret
        return ret

    def get_labels(self):
        self._write_cif()
        self._write_supercell_cif()
        self._get_mode_details()
        self.get_modes()
        return self.projection()

    def write_result(self, fname):
        self.get_labels()
        with open(fname, 'w') as myfile:
            myfile.write("# %7s  %5s  %10s  %5s  %9s\n"%("qname",'id', 'freq(cm-1)', 'label','symmetry'))
            for key, val in self.mode_info.items():
                qname, imode =key
                freq=val['freq']
                label=val['label']
                sym=val['symmetry']
                myfile.write("%8s  %5d  %10.5e  %5s  %9s\n"%(qname, imode, freq, label, sym))


def test():
    l = LabelDDB('./phbst/run.abo_PHBST.nc', tmpdir='./tmp')
    l.write_result('result.txt')


test()
