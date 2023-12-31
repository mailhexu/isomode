# Dependencies

* ASE (atomic simulation environment)
* numpy 
* spglib (for space group)
* abipy (for reading netcdf files)
* request  (for interacting with the findsym and isodistort server)
* bs4 (for parsing the html pages)
* optional:
  * anaddb compiled with netcdf (not needed by the package but it's needed for generating the *PHBST.nc file) 

## Installation

* Dependencies should be automatically installed if pip works, otherwise they need to be installed manually.  Abipy could be difficult though...

* download or git clone package 

* ```
  cd isomode
  ```

  ​

* run setup.py

  ```
  python setup.py setup --user
  ```

  or 

  ```
  python setpu.py develop --user
  ```

   (which could be useful if you want to update the package through git and don't want to reinstall every time)

  ​



## Usage

### Generating distorted structures

* prepare the phonon band structure netcdf file, in which the qpoints are the zone-center and high-symmetry points. 

Note:

​           (a) that symmetry equivalent qpoints are needed. 

​	    (b) do not add points between them. (set ndivsm to 2)

   Here's an example of the anaddb input for Dion-Jacobson structure.

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

   The (prefix)_PHBST.nc file is the one will be needed. (which has the phonon eigen vectors.)

* Make a copy of gen_all.py from the template directory or here, and modify them to your needs. The meaning of the parameters are in the comments.

  ```python
  #!/usr/bin/env python
  import numpy as np
  import os
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
          path='tmp',  # temporary directory, but perhaps you may find things useful in it?
          supercell_matrix=np.eye(3) * 2,  # supercell matrix
          max_freq=0.0,  # maximum frequency. use 0.0 if only unstable mode is required
          amp=0.03,  # amplitude of each mode
          pickle_fname='all_modes.pickle',  # output to this pickle file
          cif_dir='all_modes',  # output cif to this directory
          primitive=True  # whether to make it primitve
  )

  ```

* Run gen_all.py

  ```
  python gen_all.py
  ```

  You will get:

  - tmp/primitive.cif : the primitve cell cif file. 
  - all_modes.pickle : which contains the structure (in ase atoms format) of all the distorted structures, and the amplitudes of  modes, the irreps labels, spacegroups for each of them.
  - all_modes directory with all the cif files of the distorted structure.

### Identify symmetry adapted modes from a high symmetry structure and a low symmetry structure.

The command view_distort.py can be used for this purpose.

```
view_distort.py --help
usage: view_distort.py [-h] [-p PARENT] [-d DISTORT] [-o OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -p PARENT, --parent PARENT
                        parent cif file name
  -d DISTORT, --distort DISTORT
                        distorted cif file name
  -o OUTPUT, --output OUTPUT
                        mode details output filename
```

For example,  a cubic (cubic_std.cif) as the parent structure, and a P21/c YNiO3 structure as the distorted structure.  

```
 >>> view_distort.py -p cubic_std.cif -d P21c.cif -o mode_details.txt
 
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

The format of the output is 

​                    [qpoint]label: amplitude_in_supercell amplitude_in_primitive(parent)_cell

Note:

- Here you can see several mode should not be in P21/c structure has non-zero amplitude. (GM4-, GM5-, X1+, X5+, M5-). They are quite small and is due to numerial error. They should be ignored.

- This can only be used with internet access because it upload the files to isodistort server (http://stokes.byu.edu/iso/isodistortform.php)

- A output file is generated (use -o option to specify it), which contains the details of the mode decompositon. This file is generated by the isodistort server. The format of the file is described at http://stokes.byu.edu/iso/isodistorthelp.php

  ​





