---
title: "Symmetry Analysis"
weight: 3
---

# Symmetry Analysis

IsoMode provides tools to identify symmetry-adapted modes by comparing a high-symmetry structure with a low-symmetry structure.

## Using view_distort.py

The `view_distort.py` command is the main tool for symmetry analysis. It takes a parent structure (high symmetry) and a distorted structure (low symmetry) as inputs and analyzes the symmetry-adapted modes that relate them.

### Command Usage

```bash
view_distort.py -p PARENT -d DISTORT -o OUTPUT
```

### Command Options

- `-p`, `--parent`: Parent (high symmetry) CIF file name
- `-d`, `--distort`: Distorted (low symmetry) CIF file name
- `-o`, `--output`: Output filename for mode details

### Example Usage

Here's an example analyzing a P21/c YNiO3 structure against its cubic parent:

```bash
view_distort.py -p cubic_std.cif -d P21c.cif -o mode_details.txt
```

### Understanding the Output

The command outputs mode amplitudes and their relative contributions. Example output:

```
[1/2,1/2,1/2]R5-:  2.9567  0.7392
  [1/2,1/2,0]M2+:  2.2424  0.5606
    [0,1/2,0]X5-:  1.6549  0.4137
[1/2,1/2,1/2]R4-:  0.4241  0.1060
[1/2,1/2,1/2]R2-:  0.2839  0.0710
  [1/2,1/2,0]M3+:  0.1099  0.0275
  [1/2,1/2,0]M5+:  0.0695  0.0174
[1/2,1/2,1/2]R3-:  0.0199  0.0050
     [0,0,0]GM4-:  0.0006  0.0001
     [0,0,0]GM5-:  0.0004  0.0001
    [0,1/2,0]X1+:  0.0004  0.0001
    [0,1/2,0]X5+:  0.0003  0.0001
  [1/2,1/2,0]M5-:  0.0003  0.0001
           Total:  4.0971  1.0243
      SPACEGROUP: P2_1/c (14)
```

The output format is:
```
[qpoint]label: amplitude_in_supercell amplitude_in_primitive(parent)_cell
```

### Important Notes

1. Small Amplitudes: Some modes may show very small non-zero amplitudes (like GM4-, GM5-, X1+, X5+, M5- in the example) due to numerical errors. These should typically be ignored.

2. Internet Requirement: This tool requires internet access as it uses the ISODISTORT server (http://stokes.byu.edu/iso/isodistortform.php) for calculations.

3. Output File Details: The detailed output file (specified by `-o`) contains comprehensive mode decomposition information. The format follows the specification described at http://stokes.byu.edu/iso/isodistorthelp.php
