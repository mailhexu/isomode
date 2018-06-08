#!/usr/bin/env python3
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
