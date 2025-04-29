---
title: Installation
weight: 1
---

# Installation


## Installation Methods

### Install from PyPI (Recommended)

The easiest way to install IsoMode is through pip:

```bash
pip install isomode --user
```

### Install from Source

To install from source:
Download the source code, extract it to the isomode directory,
```bash
cd isomode
pip install . --user
```

This method is useful if you want to modify the source code or contribute to development.


### Dependencies
  IsoMode has several dependencies that need to be installed, but they are automatically installed when you install IsoMode via pip. However, if you want to install them manually, they are listed below:

## Dependencies

* ASE (atomic simulation environment)
* numpy 
* spglib (for space group)
* abipy (for reading netcdf files)
* request (for interacting with the findsym and isodistort server)
* bs4 (for parsing the html pages)

Optional dependencies:
* anaddb compiled with netcdf (not needed by the package but required for generating *PHBST.nc files)

