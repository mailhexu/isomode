#!/bin/bash

# Example shell script demonstrating the CLI usage of anamode.py

# Basic usage with default parameters (2x2x2 supercell, only Gamma point)
python -m isomode.anamode ./LAO_PHBST.nc mode_labels_basic.txt

# Advanced usage with custom parameters
# - 3x3x3 supercell
# - Multiple q-points (Gamma, X, and M points)
# - Custom temporary directory
python -m isomode.anamode \
    ./LAO_PHBST.nc \
    mode_labels_advanced.txt \
    --sc_mat 2 2 2 \
    --tmpdir ./tmp \
    --qpoints Gamma 0.0 0.0 0.0 M 0.5 0.5 0.0

echo "Mode labeling complete. Results written to mode_labels_*.txt"
