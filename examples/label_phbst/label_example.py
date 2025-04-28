#!/usr/bin/env python
import numpy as np
from isomode.anamode import label_phbst

# Define q-points to analyze
qpoints = {
    'Gamma': [0.0, 0.0, 0.0],
    #'X': [0.5, 0.0, 0.0],
    'M': [0.5, 0.5, 0.0]
}

# Define a 2x2x2 supercell matrix
sc_matrix = np.eye(3) * 2

# Example usage of label_phbst function
label_phbst(
    fname_phbst='./LAO_PHBST.nc',  # Input PHBST.nc file
    output_fname='mode_labels.txt',            # Output file for results
    sc_mat=sc_matrix,                         # Supercell matrix
    qdict=qpoints,                            # Q-points to analyze
    tmpdir='tmp'                              # Temporary directory
)

print("Mode labeling complete. Results written to mode_labels.txt")
