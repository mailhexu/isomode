# Label PHBST Examples

This directory contains examples demonstrating how to label phonon modes from ABINIT's PHBST.nc files using ISODISTORT.

## Overview

The `label_phbst` functionality helps analyze phonon modes by:
- Reading phonon data from ABINIT PHBST.nc files
- Matching modes with symmetry-adapted modes from ISODISTORT
- Providing symmetry labels and irreducible representations
- Supporting analysis at multiple q-points

## Examples

### Python API Usage (label_example.py)

The Python script demonstrates programmatic usage of the `label_phbst` function:

```python
from isomode.anamode import label_phbst
import numpy as np

# Define q-points to analyze
qpoints = {
    'Gamma': [0.0, 0.0, 0.0],
    'X': [0.5, 0.0, 0.0],
    'M': [0.5, 0.5, 0.0]
}

# Define supercell matrix (2x2x2)
sc_matrix = np.eye(3) * 2

# Label the modes
label_phbst(
    fname_phbst='path/to/run.abo_PHBST.nc',  # Input file
    output_fname='mode_labels.txt',           # Output file
    sc_mat=sc_matrix,                         # Supercell matrix
    qdict=qpoints,                            # Q-points
    tmpdir='tmp'                              # Temporary directory
)
```

### Command Line Usage (label_cli.sh)

The shell script shows two ways to use the CLI:

1. Basic usage (default parameters):
```bash
python -m isomode.anamode input.nc output.txt
```

2. Advanced usage (custom parameters):
```bash
python -m isomode.anamode \
    input.nc \
    output.txt \
    --sc_mat 2 2 2 \
    --tmpdir ./tmp \
    --qpoints \
        Gamma 0.0 0.0 0.0 \
        X 0.5 0.0 0.0 \
        M 0.5 0.5 0.0
```

## Input Parameters

### Python API Parameters

- `fname_phbst` (str): Path to ABINIT PHBST.nc file
- `output_fname` (str): Path for output results file
- `sc_mat` (numpy.ndarray): 3x3 supercell matrix (default: 2x2x2)
- `qdict` (dict): Dictionary of q-points to analyze
- `tmpdir` (str, optional): Directory for temporary files

### CLI Parameters

- `input`: Input PHBST.nc file path
- `output`: Output results file path
- `--sc_mat`: Three integers defining supercell (default: 2 2 2)
- `--tmpdir`: Temporary directory path
- `--qpoints`: Q-points as alternating labels and coordinates

## Output Format

The output file contains tabulated mode information:

```
# qname   id  freq(cm-1)  label  symmetry
Gamma     0   100.0       GM1+   A1g
Gamma     1   150.0       GM2-   Eg
X         0   200.0       X5+    B1
...
```

Each row includes:
- `qname`: Q-point label
- `id`: Mode index
- `freq`: Frequency in cm-1
- `label`: Mode label
- `symmetry`: Irreducible representation

## Key Features

- Analysis of multiple q-points simultaneously
- Custom supercell matrix support
- Temporary file management
- Integration with ISODISTORT for symmetry analysis
- Both programmatic and command-line interfaces
- Detailed mode labeling with frequencies and symmetries
