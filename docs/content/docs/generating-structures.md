---
title: "Generating Distorted Structures"
weight: 4
---

# Generating Distorted Structures

IsoMode can generate distorted structures based on phonon calculations. This feature requires proper setup of phonon band structure calculations.

## Prerequisites

You need a phonon band structure netcdf file that includes:
- Zone-center points
- High-symmetry points
- Symmetry equivalent qpoints

Important notes for phonon calculations:
- Symmetry equivalent qpoints must be included
- Do not add points between high-symmetry points (set ndivsm to 2)

### Example ANADDB Input

Here's an example input for a Dion-Jacobson structure:

```
# ANADB input for phonon bands and DOS
ndivsm 2
nqpath 8
qpath
   0.0    0.0    0.0 # Gamma
   0.0    0.5    0.0 # Xy
   0.5    0.0    0.0 # Xx
   0.5    0.5    0.0 # M
   0.5    0.0    0.5 # Rx
   0.0    0.5    0.5 # Ry
   0.5    0.5    0.5 # A
   0.0    0.0    0.5 # Z
asr 2
ngqpt 2 2 2
chneut 1
dipdip 0
ifcflag 1
nqshft 1
q1shft 0 0 0
```

## Using gen_all.py

1. Create a copy of `gen_all.py` from the template directory
2. Modify the parameters according to your needs
3. Run the script

### Example Script

```python
#!/usr/bin/env python
import numpy as np
from isomode.gen_all import run_all

if __name__=="__main__":
    run_all(
        fname='./run.abo_PHBST.nc',  # phonon band netcdf file
        qdict={
            'Gamma': [0.0, 0.0, 0.0],
            'Xy': [0, 0.5, 0],
            'Xx': [0.5, 0.0, 0],
            'M': [0.5, 0.5, 0],
            'Rx': [0.5, 0.0, 0.5],
            'Ry': [0.0, 0.5, 0.5],
            'A': [0.5, 0.5, 0.5],
            'Z': [0, 0, 0.5]
        },  # qpoints in netcdf file
        path='tmp',  # temporary directory
        supercell_matrix=np.eye(3) * 2,  # supercell matrix
        max_freq=0.0,  # maximum frequency (0.0 for only unstable modes)
        amp=0.03,  # amplitude of each mode
        pickle_fname='all_modes.pickle',  # output pickle file
        cif_dir='all_modes',  # output directory for CIF files
        primitive=True  # whether to make it primitive
    )
```

### Running the Script

To generate the structures:

```bash
python gen_all.py
```

## Output Files

The script generates several outputs:

1. `tmp/primitive.cif`: The primitive cell CIF file

2. `all_modes.pickle`: Contains:
   - Distorted structures (in ASE atoms format)
   - Mode amplitudes
   - Irreps labels
   - Spacegroups for each structure

3. `all_modes/` directory: Contains CIF files for all generated distorted structures

### Parameters Explanation

- `fname`: Path to the PHBST netcdf file
- `qdict`: Dictionary mapping q-point names to coordinates
- `path`: Directory for temporary files
- `supercell_matrix`: Matrix defining the supercell
- `max_freq`: Maximum frequency cutoff (0.0 for unstable modes only)
- `amp`: Amplitude for mode distortions
- `pickle_fname`: Output file for structural data
- `cif_dir`: Directory for output CIF files
- `primitive`: Whether to generate primitive cells
