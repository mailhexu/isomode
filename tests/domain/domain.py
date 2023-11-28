import numpy as np
import sys
# from phonopy.structure.cells import get_supercell
from ase import Atoms
from ase.io import write

import numpy as np
import copy
import pyDFTutils.perovskite.perovskite_mode as perovskite_mode
from pyDFTutils.perovskite.perovskite_mode import Gamma_modes

import spglib.spglib
#from phonopy.structure.atoms import PhonopyAtoms as Atoms
from pyDFTutils.perovskite.cubic_perovskite import gen_primitive
from pyDFTutils.ase_utils import vesta_view
from isomode.frozen_mode import distorted_cell, default_mod_func
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from functools import partial


def gen_domain(
        name,
        cell=[3.9, 3.9, 3.9],
        supercell_matrix=[[1, -1, 0], [1, 1, 0], [0, 0, 2]],
        modes=dict(
            R2_m_O1=(0.0, default_mod_func,
                     0.0),  # R2-[O1:c:dsp]A2u(a), O, breathing
            R3_m_O1=(
                0.0, default_mod_func, 0.0
            ),  # R3-[O1:c:dsp]A2u(a), O JT inplane-stagger, out-of-plane antiphase
            R3_m_O2=(
                0.0, default_mod_func, 0.0
            ),  # R3-[O1:c:dsp]A2u(b), O, out-of-plane-stagger, inplane antiphase, Unusual.
            R4_m_A1=(0.0, default_mod_func,
                     0.0),  # R4-[Nd1:a:dsp]T1u(a), A , Unusual
            R4_m_A2=(0.0, default_mod_func,
                     0.0),  # R4-[Nd1:a:dsp]T1u(b), A, Unusual
            R4_m_A3=(0.0, default_mod_func,
                     0.0),  # R4-[Nd1:a:dsp]T1u(c), A, Unusual
            R4_m_O1=(0.0, default_mod_func,
                     0.0),  # R4-[O1:c:dsp]Eu(a), O, Unusual
            R4_m_O2=(0.0, default_mod_func,
                     0.0),  # R4-[O1:c:dsp]Eu(b), O, Unusual
            R4_m_O3=(0.0, default_mod_func,
                     0.0),  # R4-[O1:c:dsp]Eu(c), O, Unusual
            R5_m_O1=(0.0, default_mod_func, 0.0),  # R5-[O1:c:dsp]Eu(a), O  a-
            R5_m_O2=(0.0, default_mod_func, 0.0),  # R5-[O1:c:dsp]Eu(b), O  b-
            R5_m_O3=(0.0, default_mod_func, 0.0),  # R5-[O1:c:dsp]Eu(c), O  c-
            #X3_m_A1=(0.0, default_mod_func,
            #         0.0),  # X3-[Nd1:a:dsp]T1u(a), What's this..
            #X3_m_O1=(0.0, default_mod_func, 0.0),  # X3-[O1:c:dsp]A2u(a)

            #X5_m_A1=(0.0, default_mod_func, 0.0),  # [Nd1:a:dsp]T1u(a), A , Antiferro mode
            #X5_m_A2=(0.0, default_mod_func, 0.0),  # [Nd1:a:dsp]T1u(b), A , save as above
            #X5_m_O1=(0.0, default_mod_func, 0.0),  # [Nd1:a:dsp]T1u(a), O , Antiferro mode
            #X5_m_O2=(0.0, default_mod_func, 0.0),  # [Nd1:a:dsp]T1u(b), O , same as above
            #M2_p_O1=(0.0, default_mod_func, 0.0),  # M2+[O1:c:dsp]Eu(a), O, In phase rotation c+
            Z5_m_A1=(0.0, default_mod_func,
                     0.0),  # [Nd1:a:dsp]T1u(a), A , Antiferro mode
            Z5_m_A2=(0.0, default_mod_func,
                     0.0),  # [Nd1:a:dsp]T1u(b), A , save as above
            Z5_m_O1=(0.0, default_mod_func,
                     0.0),  # [Nd1:a:dsp]T1u(a), O , Antiferro mode
            Z5_m_O2=(0.0, default_mod_func,
                     0.0),  # [Nd1:a:dsp]T1u(b), O , same as above
            M2_p_O1=(0.0, default_mod_func,
                     0.0),  # M2+[O1:c:dsp]Eu(a), O, In phase rotation
            M3_p_O1=(0.0, default_mod_func,
                     0.0),  # M3+[O1:c:dsp]A2u(a), O, D-type JT inplane stagger
            M5_p_O1=(0.0, default_mod_func,
                     0.0),  # M5+[O1:c:dsp]Eu(a), O, Out of phase tilting
            M5_p_O2=(0.0, default_mod_func,
                     0.0),  # M5+[O1:c:dsp]Eu(a), O, Out of phase tilting
            M4_p_O1=(
                0.0,
                default_mod_func, 0.0
            ),  # M4+[O1:c:dsp]A2u(a), O, in-plane-breathing (not in P21/c)
            G_Ax=(0.0, default_mod_func, 0.0),
            G_Ay=(0.0, default_mod_func, 0.0),
            G_Az=(0.0, default_mod_func, 0.0),
            G_Sx=(0.0, default_mod_func, 0.0),
            G_Sy=(0.0, default_mod_func, 0.0),
            G_Sz=(0.0, default_mod_func, 0.0),
            G_Axex=(0.0, default_mod_func, 0.0),
            G_Axey=(0.0, default_mod_func, 0.0),
            G_Axez=(0.0, default_mod_func, 0.0),
            G_Lx=(0.0, default_mod_func, 0.0),
            G_Ly=(0.0, default_mod_func, 0.0),
            G_Lz=(0.0, default_mod_func, 0.0),
            G_G4x=(0.0, default_mod_func, 0.0),
            G_G4y=(0.0, default_mod_func, 0.0),
            G_G4z=(0.0, default_mod_func, 0.0),
        ), return_nondistorted=True):
    atoms = gen_primitive(name=name, mag_order='PM', latticeconstant=cell[0])
    spos = atoms.get_scaled_positions()
    atoms.set_cell(cell)
    atoms.set_scaled_positions(spos)
    dcell = distorted_cell(atoms, supercell_matrix=supercell_matrix)
    eigvec = np.zeros(15)

    mode_dict = {
        'R2_m_O1': perovskite_mode.R2p,
        'R3_m_O1': perovskite_mode.
        R12p_1,  # R3-[O1:c:dsp]A2u(a), O JT inplane-stagger, out-of-plane antiphase
        'R3_m_O2': perovskite_mode.
        R12p_2,  # R3-[O1:c:dsp]A2u(b), O, out-of-plane-stagger, inplane antiphase
        'R4_m_A1': perovskite_mode.R15_1,  # R4-[Nd1:a:dsp]T1u(a), A
        'R4_m_A2': perovskite_mode.R15_2,  # R4-[Nd1:a:dsp]T1u(b), A
        'R4_m_A3': perovskite_mode.R15_3,  # R4-[Nd1:a:dsp]T1u(c), A
        'R4_m_O1': perovskite_mode.R15_4,  # R4-[O1:c:dsp]Eu(a), O
        'R4_m_O2': perovskite_mode.R15_5,  # R4-[O1:c:dsp]Eu(b), O
        'R4_m_O3': perovskite_mode.R15_6,  # R4-[O1:c:dsp]Eu(c), O
        'R5_m_O1': perovskite_mode.
        R25_1,  # R5-[O1:c:dsp]Eu(a), O, out-of-phase rotation a-
        'R5_m_O2': perovskite_mode.R25_2,  # R5-[O1:c:dsp]Eu(b), O, b-
        'R5_m_O3': perovskite_mode.
        R25_3,  # R5-[O1:c:dsp]Eu(c), O, c-. For Pnma. Do not use.

        #'X3_m_A1':perovskite_mode., # X3-[Nd1:a:dsp]T1u(a), What's this..
        #'X3_m_O1':perovskite_mode., # X3-[O1:c:dsp]A2u(a)
        'Z5_m_A1':
        perovskite_mode.Z5p_1,  # [Nd1:a:dsp]T1u(a), A , Antiferro mode
        'Z5_m_A2':
        perovskite_mode.Z5p_2,  # [Nd1:a:dsp]T1u(b), A , save as above
        'Z5_m_O1':
        perovskite_mode.Z5p_3,  # [Nd1:a:dsp]T1u(b), O , same as above
        'Z5_m_O2':
        perovskite_mode.Z5p_4,  # [Nd1:a:dsp]T1u(b), O , same as above
        'M2_p_O1':
        perovskite_mode.M3,  # M2+[O1:c:dsp]Eu(a), O, In phase rotation
        'M3_p_O1': perovskite_mode.
        M2,  # M3+[O1:c:dsp]A2u(a), O, D-type JT inplane stagger
        'M5_p_O1':
        perovskite_mode.M5_1,  # M5+[O1:c:dsp]Eu(a), O, Out of phase tilting
        'M5_p_O2': perovskite_mode.
        M5_2,  # M5+[O1:c:dsp]Eu(b), O, Out of phase tilting, -above
        'M4_p_O1': perovskite_mode.
        M4,  # M4+[O1:c:dsp]A2u(a), O, in-plane-breathing (not in P21/c)
    }

    # add Gamma modes to mode_dict
    Gamma_mode_dict = Gamma_modes(atoms.get_chemical_symbols())
    mode_dict.update(Gamma_mode_dict)

    mode_disps = {}
    qdict = {
        'G': [0, 0, 0],
        #'X':[0,0.0,0.5],
        'M': [0.5, 0.5, 0],
        'R': [0.5, 0.5, 0.5],
        'Z': [0.0, 0.0, 0.5]
    }
    disps = 0.0  #np.zeros(3,dtype='complex128')
    for key, val in modes.items():
        name = key
        amp, mod_func, argument = val
        eigvec = np.array(mode_dict[name])
        disp = dcell._get_displacements(
            eigvec=eigvec,
            q=qdict[name[0]],
            amplitude=amp,
            argument=0,
            mod_func=mod_func)
        disps += disp

    nondistorted = dcell._supercell
    nondistorted.set_pbc(True)
    newcell = dcell._get_cell_with_modulation(disps)
    newcell = Atoms(newcell)
    print(spglib.get_spacegroup(newcell))
    if return_nondistorted:
        return newcell, nondistorted
    else:
        return newcell


def isotropy_normfactor(scell, sc_mat, disps):
    """
    pcell: primitive cell parameter. 3*3 matrix
    sc_mat: primitive->supercell transform matrix. 3*3 matrix
    disps: list of vectors defining displacements.
    """
    sum = 0.0
    for disp in disps:
        sum += (np.linalg.norm(np.dot(scell, disp)))**2
        norm_factor = 1.0 / np.sqrt(sum)
    return norm_factor


def conv_circ(signal, ker):
    '''
        signal: real 1D array
        ker: real 1D array
        signal and ker must have same shape
    '''
    return np.real(np.fft.ifft(np.fft.fft(signal) * np.fft.fft(ker)))


from scipy.ndimage import convolve


def smooth_domain_func(R, width, ncell, pos, dim, d0_value=1.0, d1_value=-1.0):
    """
    This funcunction create a smoothed step function in x, y, or z direction. The size is ncell. The value is 1 in a range (0~pos*ncell), and low_value in other part.
    ncell: number of cells along one direction
    width: domain wall width in the unit of cell.
    pos: position of the domain wall . Note that there are two domains, so there is alwasys another
         domain wall at 0. 0<pos<1.
    dim: 0,1,2-> x, y, or z direction.
    d0_value: amplitude coefficient in first domain
    d1_value: amplitude coefficient in second domain
    """
    r = 100
    x0 = np.linspace(0, ncell, ncell * r, endpoint=True)
    pos = int(np.floor(ncell * r * pos))
    y0 = np.zeros_like(x0)
    y0[:pos] = d0_value
    y0[pos:] = d1_value
    y0 = convolve(y0, [1.0 / r / width] * r * width, mode='wrap')
    return interp1d(x0, y0)(R[:, dim])


def mod_func_1d(R, width, ncell, pos, dim):
    return smooth_domain_func(width, ncell)(R[:, dim])


mod_func_x = partial(smooth_domain_func, width=4, ncell=16, pos=0.5, dim=0)
mod_func_y = partial(smooth_domain_func, width=4, ncell=16, pos=0.5, dim=1)
mod_func_z = partial(smooth_domain_func, width=4, ncell=16, pos=0.5, dim=2)

mod_func_x10 = partial(
    smooth_domain_func, width=4, ncell=16, pos=0.5, dim=0, d1_value=0.0)
mod_func_y10 = partial(
    smooth_domain_func, width=4, ncell=16, pos=0.5, dim=1, d1_value=0.0)
mod_func_z10 = partial(
    smooth_domain_func, width=4, ncell=16, pos=0.5, dim=2, d1_value=0.0)

mod_func_x01 = partial(
    smooth_domain_func,
    width=4,
    ncell=16,
    pos=0.5,
    dim=0,
    d0_value=0.0,
    d1_value=1.0)
mod_func_y01 = partial(
    smooth_domain_func,
    width=4,
    ncell=16,
    pos=0.5,
    dim=1,
    d0_value=0.0,
    d1_value=1.0)
mod_func_z01 = partial(
    smooth_domain_func,
    width=4,
    ncell=16,
    pos=0.5,
    dim=2,
    d0_value=0.0,
    d1_value=1.0)


def antipolar_180_x():
    atoms, orig = gen_domain(
        name='PbTiO3',
        cell=[3.7, 3.7, 3.7],
        supercell_matrix=[[16, 0, 0], [0, 2, 0], [0, 0, 2]],
        modes=dict(
            Z5_m_A1=(1.5, mod_func_x,
                     0.0),  # [Nd1:a:dsp]T1u(a), A , Antiferro mode
            Z5_m_A2=(1.5, mod_func_x,
                     0.0),  # [Nd1:a:dsp]T1u(b), A , save as above
            Z5_m_O1=(2.0, mod_func_x,
                     0.0),  # [Nd1:a:dsp]T1u(a), O , Antiferro mode
            Z5_m_O2=(2.0, mod_func_x,
                     0.0))  # [Nd1:a:dsp]T1u(b), O , same as above

    )
    atoms.set_pbc([True, True, True])
    write('antipolar_180domain_x_nondistorted.cif', orig)
    write('antipolar_180domain_x.cif', atoms)
    vesta_view(atoms)


def antipolar_180_z():
    # amplitude, mod_func
    atoms, orig = gen_domain(
        name='PbTiO3',
        cell=[3.7, 3.7, 3.7],
        supercell_matrix=[[2, 0, 0], [0, 2, 0], [0, 0, 16]],
        modes=dict(
            R5_m_O1=(2.0, default_mod_func, 0.0),  # R5-[O1:c:dsp]Eu(a), O  a-
            R5_m_O2=(2.0, default_mod_func, 0.0),  # R5-[O1:c:dsp]Eu(b), O  b-
            R5_m_O3=(2.0, default_mod_func, 0.0),  # R5-[O1:c:dsp]Eu(c), O  c-
 
            Z5_m_A1=(1.5, mod_func_z,
                     0.0),  # [Nd1:a:dsp]T1u(a), A , Antiferro mode
            Z5_m_A2=(1.5, mod_func_z,
                     0.0),  # [Nd1:a:dsp]T1u(b), A , save as above
            Z5_m_O1=(2.0, mod_func_z,
                     0.0),  # [Nd1:a:dsp]T1u(a), O , Antiferro mode
            Z5_m_O2=(2.0, mod_func_z, 0.0),

        ))  
    atoms.set_pbc([True, True, True])

    write('antipolar_180domain_z_nondistorted.cif', orig)
    write('antipolar_180domain_z.cif', atoms)
    vesta_view(atoms)


def polar_180_x_longi():
    # amplitude, mod_func
    atoms, orig = gen_domain(
        name='PbTiO3',
        cell=[3.901, 4.128, 3.901],
        supercell_matrix=[[16, 0, 0], [0, 1, 0], [0, 0, 1]],
        modes=dict(
            G_Ax=(0.0, default_mod_func, 0.0),
            G_Ay=(0.0, default_mod_func, 0.0),
            G_Az=(0.0, default_mod_func, 0.0),
            G_Sx=(0.0, default_mod_func, 0.0),
            G_Sy=(1.0, mod_func_x, 0.0),
            G_Sz=(0.0, mod_func_x, 0.0),
            G_Axex=(0.0, default_mod_func, 0.0),
            G_Axey=(0.0, default_mod_func, 0.0),
            G_Axez=(0.0, default_mod_func, 0.0),
            G_Lx=(0.0, default_mod_func, 0.0),
            G_Ly=(2.0, mod_func_x, 0.0),
            G_Lz=(0.0, mod_func_x, 0.0),
            G_G4x=(0.0, default_mod_func, 0.0),
            G_G4y=(0.0, default_mod_func, 0.0),
            G_G4z=(0.0, default_mod_func, 0.0),
        ))
    atoms.set_pbc([True, True, True])
    write('polar_180domain_nondistored.cif', orig)
    write('polar_180domain.cif', atoms)
    vesta_view(atoms)


def polar_180_x_trans():
    atoms, orig = gen_domain(
        name='PbTiO3',
        cell=[3.901, 3.7, 3.7],
        supercell_matrix=[[16, 0, 0], [0, 2, 0], [0, 0, 2]],
        modes=dict(
            G_Sx=(3.0, mod_func_x, 0.0),
        ))
    atoms.set_pbc(True)
    write("polar_180_x_nondistorted.cif",atoms)
    write("polar_180_x.cif",atoms)
    vesta_view(atoms)

def polar_180_x_trans_phonon():
    phon = PhonMode.loadDDB('PbTiO3.ddb')
    atoms= phon.gen_domain(
        supercell_matrix=[[16, 0, 0], [0, 2, 0], [0, 0, 2]],
        modes={((0.0,0.0,0.0), 1): (3.0, mod_func_x, 0.0)}
        )
    atoms.set_pbc(True)
    write("polar_180_x_phonon.cif",atoms)
    vesta_view(atoms)


def polar_90_xy_trans():
    # amplitude, mod_func
    atoms, orig = gen_domain(
        name='PbTiO3',
        cell=[3.7, 3.7, 3.7],
        supercell_matrix=[[16, 0, 0], [0, 2, 0], [0, 0, 2]],
        modes=dict(
            G_Ax=(0.0, default_mod_func, 0.0),
            G_Ay=(0.0, default_mod_func, 0.0),
            G_Az=(0.0, default_mod_func, 0.0),
            G_Sx=(3.0, mod_func_x10, 0.0),
            G_Sy=(3.0, mod_func_x01, 0.0),
            G_Sz=(0.0, default_mod_func, 0.0),
            G_Axex=(0.0, default_mod_func, 0.0),
            G_Axey=(0.0, default_mod_func, 0.0),
            G_Axez=(0.0, default_mod_func, 0.0),
            G_Lx=(0.0, default_mod_func, 0.0),
            G_Ly=(0.0, default_mod_func, 0.0),
            G_Lz=(0.0, default_mod_func, 0.0),
            G_G4x=(0.0, default_mod_func, 0.0),
            G_G4y=(0.0, default_mod_func, 0.0),
            G_G4z=(0.0, default_mod_func, 0.0),
        ))
    atoms.set_pbc(True)
    write("polar_90_nondistorted.cif", orig)
    write("polar_90.cif", atoms)
    vesta_view(atoms)

def cyc_func_xy(R):
    return np.sin(R[:, 1]*4*np.pi/16.0) #+ np.cos(R[:, 1]*4*np.pi/16.0)


def cyc_func_yx(R):
    return np.cos(R[:, 0]*4*np.pi/16.0) #+ np.sin(R[:, 1]*4*np.pi/16.0)


def polar_cycle_xy():
    # amplitude, mod_func
    atoms, orig = gen_domain(
        name='PbTiO3',
        cell=[3.7, 3.7, 3.7],
        supercell_matrix=[[16, 0, 0], [0, 16, 0], [0, 0, 2]],
        modes=dict(
            G_Ax=(0.0, default_mod_func, 0.0),
            G_Ay=(0.0, default_mod_func, 0.0),
            G_Az=(0.0, default_mod_func, 0.0),
            G_Sx=(9.0, cyc_func_xy, 0.0),
            G_Sy=(9.0, cyc_func_yx, 0.0),
            G_Sz=(0.0, default_mod_func, 0.0),
            G_Axex=(0.0, default_mod_func, 0.0),
            G_Axey=(0.0, default_mod_func, 0.0),
            G_Axez=(0.0, default_mod_func, 0.0),
            G_Lx=(0.0, default_mod_func, 0.0),
            G_Ly=(0.0, default_mod_func, 0.0),
            G_Lz=(0.0, default_mod_func, 0.0),
            G_G4x=(0.0, default_mod_func, 0.0),
            G_G4y=(0.0, default_mod_func, 0.0),
            G_G4z=(0.0, default_mod_func, 0.0),
        ))
    write("polar_cycle_nondistorted.cif", orig)
    write("polar_cycle.cif", atoms)
    vesta_view(atoms)



antipolar_180_x()
antipolar_180_z()
polar_180_x_longi()
polar_90_xy_trans()
polar_cycle_xy()
