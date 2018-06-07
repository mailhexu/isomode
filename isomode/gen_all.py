#!/usr/bin/env python3
import numpy as np
import os
from isomode.phonmode import phonon_distort_generator
from isomode.isodistort_parser import Isomode
from isomode.pydistort import isodistort, view_distort


def run_all(
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
):
    dg = phonon_distort_generator(fname, qdict)

    # generate primitive cell struture. write to path/primitive.cif
    dg.generate_primitive_structure(path=path)

    # generating single phonon mode structure
    # comment if not needed
    ds = dg.generate_distorted_structures(
        supercell_matrix=supercell_matrix,
        amplitude=0.5,
        unstable_only=True,
        path=os.path.join(path, 'single_modes'))

    # generate structure with all (unstable) phonons.
    # structure is written to
    ds = dg.generate_combined_mode_distorted_structures(
        supercell_matrix=supercell_matrix, max_freq=max_freq, path=path)

    # use isodistort to analyze the structure with all phonons.
    #iso = isodistort(
    #    parent_cif=os.path.join(path, 'primitive.cif'),
    #    distorted_cif=os.path.join(path, 'combined_modes.cif'))
    #ampt = iso.get_mode_amplitude_text()
    #iso.get_mode_amplitude_text()

    iso = view_distort(
        parent_fname=os.path.join(path, 'primitive.cif'),
        distorted_fname=os.path.join(path, 'combined_modes.cif'),
        out_fname=os.path.join(path, 'mode_detail.txt'))

    det_txt = os.path.join(path,'mode_detail.txt')
    #det_txt='mode_detail.txt'
    #mode_details = iso.get_mode_details(save_fname=det_txt)

    # read mode definition and amplitudes from 'mode_detail.txt'
    myparser = Isomode(det_txt)
    myparser.read_mode_amplitudes(det_txt)
    myparser.prepare_structure(
        pickle_fname=pickle_fname,  # file name of pickle file
        cif_dir=cif_dir,  # cif file directory
        amp=amp,  # amplitude
        primitive=primitive# whether to make it primitve
    )

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
