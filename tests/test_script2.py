"""Test the improved distortion analysis"""
import os
import tempfile
from isomode.pydistort import isodistort, isocif
from ase import Atoms
import numpy as np
import logging
import time

# Setup logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

# Create test directory
TEST_DIR = os.getcwd() #tempfile.mkdtemp()
print(f"Using test directory: {TEST_DIR}")

# Create a simple test structure - cubic perovskite
# SrTiO3 structure: Sr at corners, Ti at center, O at face centers
a = 3.9  # lattice parameter in Å
atoms = Atoms('SrTiO3',
              positions=[[0, 0, 0],         # Sr at origin
                        [a/2, a/2, a/2],    # Ti at center
                        [a/2, a/2, 0],      # O at face centers
                        [a/2, 0, a/2],
                        [0, a/2, a/2]],
              cell=[a, a, a],
              pbc=True)

# First get high symmetry structure
print("\nGetting high symmetry structure...")
parent_cif = os.path.join(TEST_DIR, 'parent.cif')
atoms.write(parent_cif)

iso_finder = isocif(parent_cif)
iso_finder.upload_cif()
iso_finder.findsym()
parent_sym_cif = os.path.join(TEST_DIR, 'parent_sym.cif')
iso_finder.save_cif(fname=parent_sym_cif)


# Create distorted structure
print("\nCreating distorted structure...")
distorted = atoms.copy()
distorted.positions[2,0] += 0.1  # Move O atom along x
distorted_cif = os.path.join(TEST_DIR, 'distorted.cif')
distorted.write(distorted_cif)

print("\nStarting distortion analysis from high-symmetry parent...")
iso = isodistort(parent_cif=parent_sym_cif, distorted_cif=distorted_cif)

# After each step, print diagnostics
#print("\nParent upload response:")
#print(iso.upload_parent_cif_text[:500])
#time.sleep(2)

#print("\nDistorted upload response:")
#print(iso.upload_distorted_cif_text[:500])
#time.sleep(2)

#print("\nBasis selection response:")
#print(iso.select_basis_text)
#time.sleep(2)

mode_details = iso.get_mode_details()

#if hasattr(iso, 'mode_details_text'):
    #print("\nMode details text:")
    #print(iso.mode_details_text[:500])

print("\nExtracted mode details:")
print(mode_details)
