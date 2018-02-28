import os
import numpy as np
from isomode.phonmode import phonon_distort_generator
from isomode.pydistort import isodistort


def test_phbst_modes(fname='./phonon/run.abo_PHBST.nc',qdict = {
        'Gamma': [0.0, 0.0, 0.0],
        'Xy': [0, 0.5, 0],
        'Xx': [0.5, 0.0, 0],
        'M': [0.5, 0.5, 0],
        'Rx': [0.5, 0.0, 0.5],
        'Ry': [0.0, 0.5, 0.5],
        'A': [0.5, 0.5, 0.5],
        'Z': [0, 0, 0.5]
    }, path='nmodes'):
    dg = phonon_distort_generator(fname, qdict)
    dg.generate_primitive_structure(path='nmodes')
    ds = dg.generate_distorted_structures(
        supercell_matrix=np.eye(3) * 2,
        amplitude=1.0,
        unstable_only=True,
        path=path)
    return ds


def gen_mode_details(parent_cif='phonon/prim_sym.cif',
                   distorted_cif='nmodes/A_0.cif',
                   mode_detail_file='mode_detail.txt'):
    iso = isodistort(parent_cif=parent_cif, distorted_cif=distorted_cif)
    ampt = iso.get_mode_amplitude_text()
    iso.get_mode_amplitude_text()
    mode_details = iso.get_mode_details(save_fname=mode_detail_file)

def gen_all():
    path='nmodes'
    ds=test_phbst_modes(path=path)
    for d in ds:
        modename='%s_%s'%(d['qname'],d['ibranch'])
        cifname=os.path.join(path, modename+'.cif')
        gen_mode_details(parent_cif='prim_sym.cif', distorted_cif=cifname, mode_detail_file=os.path.join('mode_details', modename+'.txt'))

gen_all()
