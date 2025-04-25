"""Test ISODISTORT mode detection"""
import os
import tempfile
from isomode.pydistort import isodistort, test_isocif
from ase import Atoms
import numpy as np

# Create test directory
TEST_DIR = tempfile.mkdtemp()
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

parent_cif = os.path.join(TEST_DIR, 'parent.cif')
atoms.write(parent_cif)

# First convert parent structure to highest symmetry
print("\nConverting to high-symmetry structure...")
test_isocif(parent_cif)
parent_sym_cif = 'save.cif'

# Create distorted structure by moving an O atom
print("\nCreating distorted structure...")
distorted = atoms.copy()
distorted.positions[2,0] += 0.1  # Move O atom along x
distorted_cif = os.path.join(TEST_DIR, 'distorted.cif')
distorted.write(distorted_cif)

# Analyze distortion
print("\nAnalyzing distortion...")
iso = isodistort(parent_cif=parent_sym_cif, distorted_cif=distorted_cif)

# Try to print complete response
print("\nDebug: Parent upload response:")
print(iso.upload_parent_cif_text[:200] + "...")

print("\nDebug: Distorted upload response:")
print(iso.upload_distorted_cif_text[:200] + "...")

print("\nDebug: Select basis response:")
print(iso.select_basis_text[:200] + "...")

print("\nDebug: Init distortion response:")
print(iso.init_distortion_text[:200] + "...")

mode_details = iso.get_mode_details()
print("\nMode details:")
print(mode_details)

mode_amps = iso.get_mode_amplitude_text()
print("\nMode amplitudes:")
print(mode_amps)
