# Anamode Documentation

This document describes the usage of the anamode module for analyzing phonon modes.

## Python Interface

### process_phbst

```python
def process_phbst(fname, tmpdir=None, sc_mat=np.eye(3) * 2, qdict=None, output_fname='result.txt')
```

Process a PHBST.nc file to analyze phonon modes and write the results to a file.

#### Parameters

- `fname` (str): Path to the PHBST.nc file
- `tmpdir` (str, optional): Temporary directory path for intermediate files
- `sc_mat` (numpy.ndarray, optional): Supercell matrix. Default is 2×2×2 diagonal matrix
- `qdict` (dict, optional): Dictionary mapping q-point names to coordinates
- `output_fname` (str, optional): Output file path. Default is 'result.txt'

#### Returns

LabelDDB: The initialized LabelDDB instance for further processing if needed

#### Examples

```python
from isomode.anamode import process_phbst

# Basic usage
process_phbst('run.abo_PHBST.nc', output_fname='output.txt')

# With custom q-point
process_phbst(
    'run.abo_PHBST.nc',
    qdict={'M': [0.5, 0.5, 0]},
    output_fname='output.txt'
)

# With custom supercell and temporary directory
process_phbst(
    'run.abo_PHBST.nc',
    tmpdir='./tmp',
    sc_mat=np.diag([3, 3, 3]),
    output_fname='output.txt'
)
```

## Command Line Interface

The anamode module can be run from the command line with various options.

### Basic Usage

```bash
python anamode.py [-h] [-c CONFIG] [-f FNAME] [-o OUTPUT] [-t TMPDIR] 
                 [-s SUPERCELL] [-q QPOINTS] [-p QPOINT]
```

### Arguments

- `-h, --help`: Show help message and exit
- `-c, --config`: TOML configuration file path
- `-f, --fname`: Input PHBST.nc file path
- `-o, --output`: Output file path
- `-t, --tmpdir`: Temporary directory path
- `-s, --supercell`: Supercell matrix (format: x,y,z, default: 2,2,2)
- `-q, --qpoints`: JSON string for q-points dictionary
- `-p, --qpoint`: Single q-point (format: name,x,y,z)

### Examples

1. Basic usage:
```bash
python anamode.py -f run.abo_PHBST.nc -o output.txt
```

2. With single q-point:
```bash
python anamode.py -f run.abo_PHBST.nc -o output.txt -p "M,0.5,0.5,0"
```

3. With custom supercell:
```bash
python anamode.py -f run.abo_PHBST.nc -o output.txt -s "3,3,3"
```

4. Using TOML configuration:
```bash
python anamode.py -c config.toml
```

## TOML Configuration

The TOML configuration file allows specifying all parameters in a structured format.

### Example Configuration

```toml
[phonon]
# Required parameters
fname = "run.abo_PHBST.nc"
output = "output.txt"

# Optional parameters
tmpdir = "./tmp"
supercell = [2, 2, 2]  # Will be converted to diagonal matrix

qpoints = { Gamma = [0.0, 0.0, 0.0], M = [0.5, 0.5, 0.0], R = [0.5, 0.5, 0.5] }
```

### Notes

- Command line arguments override values from the config file
- When using a config file, fname and output can be specified in either the config file or command line
- The supercell parameter in the TOML file is automatically converted to a diagonal matrix