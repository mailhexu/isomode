import numpy as np
from ase import Atoms
from ase.io import read, write
import os
import tempfile
try:
    from isomode.pydistort import isocif, isodistort
except ImportError:
    print("Warning: pydistort not available")

class LabelPhbst():
    def __init__(self, fname, tmpdir=None):
        self.fname = fname
        self.mode_fname = None
        self.mode_info = {}
        self.parent_atoms = None
        self.qpoints = None
        self.freqs = None
        self.displacements = None
        self.qdict = {}
        self.efreqs = {}
        self.edisps = {}
        self.idq = {}
        if tmpdir:
            self.tmpdir = tmpdir
        self.read_phonon()

    def read_phonon(self):
        """Read phonon data"""
        from abipy.dfpt.phonons import PhononBands
        pb = PhononBands.from_file(self.fname)
        struct = pb.structure

        # Convert abipy Structure to ASE Atoms
        symbols = [specie.symbol for specie in struct.species]
        positions = struct.frac_coords
        cell = struct.lattice.matrix

        self.parent_atoms = Atoms(
            symbols=symbols,
            scaled_positions=positions,
            cell=cell,
            pbc=True
        )
        
        self.qpoints = pb.qpoints.frac_coords
        self.freqs = pb.phfreqs
        self.displacements = pb.phdispl_cart
        for i, qpt in enumerate(self.qpoints):
            qpt = tuple(qpt)
            if qpt in [(0, 0, 0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.0)]:
                idir = len(self.qdict)
                if qpt == (0, 0, 0):
                    self.qdict['G'] = qpt
                    self.idq['G'] = i
                elif qpt == (0.5, 0.5, 0.0):
                    self.qdict['M'] = qpt
                    self.idq['M'] = i
                elif qpt == (0.5, 0.0, 0.0):
                    self.qdict['X'] = qpt
                    self.idq['X'] = i
                self.efreqs[idir] = self.freqs[i]
                self.edisps[idir] = self.displacements[i]

    def _write_cif(self):
        """Write structure to CIF file"""
        if not hasattr(self, 'tmpdir'):
            self.tmpdir = tempfile.mkdtemp()
        self.parent_cif = os.path.join(self.tmpdir, 'parent.cif')
        write(self.parent_cif, self.parent_atoms)

    def _get_mode_details(self):
        """Get mode details from ISODISTORT """
        if not hasattr(self, 'parent_cif'):
            self._write_cif()
        isosym = isocif(self.parent_cif)
        isosym.upload_cif()
        isosym.findsym()
        self.parent_sym_cif = os.path.join(self.tmpdir, 'parent_sym.cif')
        iso = isodistort(parent_cif=self.parent_sym_cif)
        text = iso.get_mode_details()
        return text

    def get_modes(self):
        """Parse mode information"""
        text = self._get_mode_details()
        mode_info = {}
        for line in text.split('\n'):
            if 'normfactor' in line:
                parts = line.split()
                label = parts[0].split('[')[2].split(']')[0]
                normfactor = float(parts[-1])
                mode_info[label] = {'normfactor': normfactor}
        return mode_info

def test():
    """Test the module"""
    l = LabelPhbst('../tests/YAlO3/phbst/run.abo_PHBST.nc')
    print("\nStructure Information:")
    print(f"Number of atoms: {len(l.parent_atoms)}")
    print(f"Lattice parameters: {l.parent_atoms.cell.cellpar()}")
    print(f"\nPhonon Information:")
    print(f"Number of q-points: {len(l.qpoints)}")
    for name, qpt in l.qdict.items():
        if name in l.idq:
            print(f"\nQ-point {name} {qpt}:")
            freqs = l.efreqs[name]
            print(f"Frequencies (cm-1): {freqs}")

if __name__ == "__main__":
    test()
