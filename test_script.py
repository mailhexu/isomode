import traceback
from isomode.anamode import LabelDDB
import numpy as np
import os
import tempfile

try:
    print("Testing anamode.py phonon reading and mode labeling functionality...")
    
    # Initialize with YAlO3 test data and mock service
    l = LabelDDB('tests/YAlO3/phbst/run.abo_PHBST.nc', use_mock=True)
    
    # Test structure reading
    print("\nStructure Information:")
    print(f"Number of atoms: {len(l.parent_atoms)}")
    print(f"Lattice parameters: {l.parent_atoms.cell.cellpar()}")
    print(f"Atomic positions:")
    for i, pos in enumerate(l.parent_atoms.get_scaled_positions()):
        print(f"Atom {i+1}: {pos}")
    
    # Print phonon frequencies
    print("\nPhonon Information:")
    print(f"Number of q-points: {len(l.qpoints)}")
    print("\nFrequencies at Gamma point:")
    qpt = [0.0, 0.0, 0.0]
    for i, q in enumerate(l.qpoints):
        if np.allclose(qpt, q):
            freqs = l.freqs[i]
            for j, freq in enumerate(freqs):
                print(f"Mode {j+1}: {freq:8.2f} cm-1")
            break
    
    # Test mode labeling using mock service
    print("\nMode labels from mock ISODISTORT service:")
    mode_details = l._get_mode_details()  # This uses mock service
    print(mode_details)
    
    # Parse and display mode information
    modes = l.get_modes()
    print("\nParsed mode information:")
    for label, info in modes.items():
        print(f"Mode {label}: normfactor = {info['normfactor']}")
    
    print('\nTest completed successfully')
except Exception as e:
    print('Error occurred:')
    print(traceback.format_exc())
