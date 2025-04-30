---
title: "Labeling Phonon Modes"
weight: 2
---

# Labeling Phonon Modes

The `label_phbst` function and CLI allow you to label phonon modes from ABINIT PHBST calculations.
The function can label the modes for the selected q-points. Note that the qpoints should be included in the PHBST file, and
the supercell matrix should be compatible with the qpoints.


## Python Function Usage

### Basic Usage

For basic usage with default parameters:

```python
from isomode.anamode import label_phbst

label_phbst('run.abo_PHBST.nc', 'results.txt')
```

### Advanced Usage

For custom supercell and q-points:

```python
import numpy as np
from isomode.anamode import label_phbst

# Define a 2x2x2 supercell
sc = np.eye(3) * 2

# Specify q-points of interest
qpts = {'Gamma': [0,0,0], 'X': [0.5,0,0]}

# Run the analysis
label_phbst('phonons.nc', 'modes.txt', sc_mat=sc, qdict=qpts)
```

## Command Line Interface

The module provides a convenient CLI for easy usage:

### Basic Usage

```bash
python -m isomode.anamode input.nc output.txt
```

### Advanced Usage

```bash
python -m isomode.anamode input.nc output.txt \
    --sc_mat 2 2 2 \
    --qpoints Gamma 0 0 0 X 0.5 0 0
```

### CLI Options

- `--tmpdir`: Temporary directory (default: auto-created)
- `--sc_mat`: Supercell matrix as 3 integers (default: 2 2 2)
  - Note: Should be compatible with the qpoints
- `--qpoints`: Q-points as alternating labels and coordinates
  - Example: "Gamma 0 0 0 X 0.5 0 0"

## Output Format

The output file follows this format:

```
# qname   id  freq(cm-1)  label  symmetry
Gamma     0   100.0       GM1+   A1g
Gamma     1   150.0       GM2-   Eg
...
```

Each line contains:
- qname: Name of the q-point
- id: Mode identifier
- freq: Frequency in cm⁻¹
- label: Mode label
- symmetry: Symmetry classification
