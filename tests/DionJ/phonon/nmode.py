import numpy as np
from scipy.optimize import minimize_scalar
from abipy.abilab import abiopen
from pyDFTutils.perovskite.frozen_mode import distorted_cell
from ase import Atoms
from ase.io import write
import spglib.spglib
import matplotlib.pyplot as plt
from functools import partial
from numpy import array


def displacement_cart_to_evec(displ_cart,
                              masses,
                              scaled_positions,
                              qpoint=None,
                              add_phase=True):
    """
    displ_cart: cartisien displacement. (atom1_x, atom1_y, atom1_z, atom2_x, ...)
    masses: masses of atoms.
    scaled_postions: scaled postions of atoms.
    qpoint: if phase needs to be added, qpoint must be given.
    add_phase: whether to add phase to the eigenvectors.
    """
    if add_phase and qpoint is None:
        raise ValueError('qpoint must be given if adding phase is needed')
    m = np.sqrt(np.kron(masses, [1, 1, 1]))
    evec = displ_cart * m
    if add_phase:
        phase = [
            np.exp(-2j * np.pi * np.dot(pos, qpoint))
            for pos in scaled_positions
        ]
        phase = np.kron(phase, [1, 1, 1])
        evec *= phase
        evec /= np.linalg.norm(evec)
    return evec


def rotate_evecs(evecs, theta):
    length = len(evecs)
    if length==3:
        print(length)
    if length != 2:
        return evecs
        #raise NotImplementedError("Only multiplicity=2 rotation implemented. Here multiplicity=%s"%length)
    #R=np.array([[np.cos(theta), -np.sin(theta) ],
    #            [np.sin(theta), np.cos(theta)]])
    #newevec=np.copy(evec)
    #for i in range(length):
    #    newevec[0+i*3:2+i*3]=np.dot( evec[0+i*3:2+i*3],R)
    u = np.array(evecs[0])
    v = np.array(evecs[1])
    newvec1 = np.cos(theta) * u - np.sin(theta) * v
    newvec2 = np.sin(theta) * u + np.cos(theta) * v
    return [newvec1, newvec2]


def project_to_x(vec):
    xs = vec[::3]
    ys = vec[1::3]
    return np.linalg.norm(xs, ord=1) - np.linalg.norm(ys, ord=1)


def test_align(a):
    thetas = np.arange(0, np.pi, 0.01)
    rxs = []
    for theta in thetas:
        #print(theta)
        rx = rotate_evecs(a, theta)[0]
        #print(rx)
        #print(project_to_x(rx))
        rxs.append(project_to_x(rx))
    plt.plot(thetas, rxs)
    plt.show()


rfunc = lambda theta, modes: -project_to_x(rotate_evecs(modes, theta)[0])


def align_degenerate_modes(modes):
    """
    modes are degenerate
    """
    #modes=np.real(modes)
    nmodes = len(modes)
    func = partial(rfunc, modes=modes)
    res = minimize_scalar(func, tol=1e-19)  #, bounds=[0, np.pi])
    theta=res.x
    new_modes=rotate_evecs(modes, theta)
    print("theta: ", theta / np.pi)
    for i in range(len(modes)):
        mi = modes[i] / np.linalg.norm(modes[i])
        nmi = new_modes[i] / np.linalg.norm(new_modes[i])
        #print("Aligning=========\n", np.real(mi).reshape((16,3)), np.real(nmi).reshape((16,3)))
    return new_modes


#
#a= [array([  2.53101288e-02,  -3.82088317e-02,  -4.86461363e-17]), array([  2.46900913e-01,  -3.72728069e-01,  -6.32567296e-18]), array([  2.46900913e-01,  -3.72728069e-01,  -3.08602526e-17]), array([ -7.22928661e-14,  -4.24936302e-02,   1.09275770e-16]), array([  2.13773443e-13,   1.25572159e-01,   1.92532429e-16]), array([  2.13905282e-13,   1.25572159e-01,   2.64991722e-16]), array([  1.45509166e-33,  -3.69736624e-33,  -1.68011693e-32]), array([ -5.86725440e-33,   2.72840948e-33,   2.03582253e-32]), array([ -6.60305144e-13,  -3.87737190e-01,   4.33930585e-17]), array([ -6.59999833e-13,  -3.87737190e-01,   4.13339061e-17]), array([ -4.71852214e-33,  -1.41947245e-32,  -8.36661187e-18]), array([ -1.55464941e-32,   9.19108395e-33,  -1.55779057e-17]), array([  2.25675740e-32,   4.41022949e-32,   8.36661187e-18]), array([  5.47739548e-32,  -1.94551410e-32,   1.55779057e-17]), array([  3.75366405e-13,   2.20589524e-01,   7.57601251e-17]), array([  3.75727227e-13,   2.20589524e-01,   1.29228876e-16])]

a = [
    array([-1.87833463e-18, -3.03705925e-19, 7.35475461e-33]),
    array([-1.29429689e-01, -1.76264658e-02, -8.63292537e-17]),
    array([1.29429689e-01, 1.76264658e-02, -9.54443562e-18]),
    array([8.56027336e-18, 1.28764580e-16, 1.19290208e-16]),
    array([1.25554614e-03, 9.21937203e-03, -3.31492986e-17]),
    array([-1.25554614e-03, -9.21937203e-03, 5.42249759e-17]),
    array([-4.28399726e-33, 8.79432863e-34, 3.40156914e-17]),
    array([6.91766327e-34, 1.94798870e-32, 7.25285576e-18]),
    array([-6.08154119e-02, -4.46562568e-01, 7.35384259e-17]),
    array([6.08154119e-02, 4.46562568e-01, -1.10793322e-17]),
    array([-1.83865500e-16, 5.59979131e-17, 1.78471407e-01]),
    array([-1.09030103e-16, 3.16500949e-16, 2.43052440e-02]),
    array([8.37074694e-16, -1.80417162e-17, -1.78471407e-01]),
    array([4.34089995e-17, -4.88475679e-16, -2.43052440e-02]),
    array([5.01347561e-03, 3.68135391e-02, -3.35885475e-17]),
    array([-5.01347561e-03, -3.68135391e-02, 1.25226772e-17])
]
a=[[1,1,0],
   [1,-1,0]]
#a = np.array(a).flatten()
print(a)
print(rotate_evecs(a, theta=np.pi/4))
print(align_degenerate_modes(a))
test_align(a)

#a=np.array([1,1,0, 1,-1,0],dtype=float)
#test_align()
#print(rotate_evec(a, np.pi/4))
#print(project_to_x(rotate_evec(a, np.pi/4)))
#print(rotate_evec(a, 0.774))
#print(a.reshape(16, 3))
#print(align_one_mode(a))


def align_all_modes(evals, evecs, tol=1e-7):
    """
    Here we assum evals are already sorted.
    """
    multi_modes = []
    multi_modes_inds = []
    new_evals = []
    new_evecs = []
    for i, evec in enumerate(evecs):
        # each time evec[i]!=evec[i-1], deal with mutli_modes and empty it,
        # then push eve[i] into multi_modes
        if np.abs(evals[i - 1] - evals[i]) > tol:
            if len(multi_modes) == 1:  # multi=1, no need to align
                new_evecs.append(multi_modes[0])
            elif len(multi_modes) > 1:  # multi, align and save
                for e in align_degenerate_modes(multi_modes):
                    new_evecs.append(e)
                if len(multi_modes) ==3:  # multi, align and save 
                    print(evals[i-1], evec)
            # clean multi_modes and save new
            multi_modes = []
            multi_modes.append(evec)
        else:
            multi_modes.append(evec)
    return new_evecs


def read_modes(fname='./run.abo_PHBST.nc'):
    qdict = {
        'Gamma': [0.0, 0.0, 0.0],
        'Xy': [0, 0.5, 0],
        'Xx': [0.5, 0.0, 0],
        'M': [0.5, 0.5, 0],
        'Rx': [0.5, 0.0, 0.5],
        'Ry': [0.0, 0.5, 0.5],
        'A': [0.5, 0.5, 0.5],
        'Z': [0, 0, 0.5]
    }
    phbst = abiopen(fname)
    atoms = phbst.structure.to_ase_atoms()
    write('primitive.cif',atoms)
    nbranch = 3 * len(atoms)
    masses = atoms.get_masses()
    scaled_positions = atoms.get_scaled_positions()
    myfile = open('unstable_modes.txt', 'w')
    for qname, qpt in qdict.items():
        #for i in range(nbranch):
        displ_carts = [
            phbst.get_phmode(qpt, i).displ_cart for i in range(nbranch)
        ]
        freqs = [phbst.get_phmode(qpt, i).freq for i in range(nbranch)]
        align_evecs = False
        align_disp = True
        if align_evecs:
            evecs = [
                displacement_cart_to_evec(
                    displ_cart,
                    masses=masses,
                    scaled_positions=scaled_positions,
                    qpoint=qpt,
                    add_phase=True) for displ_cart in displ_carts
            ]

            nevecs = align_all_modes(freqs, evecs)
        elif align_disp:
            ndispl_carts = align_all_modes(freqs, displ_carts)
            nevecs = [
                displacement_cart_to_evec(
                    displ_cart,
                    masses=masses,
                    scaled_positions=scaled_positions,
                    qpoint=qpt,
                    add_phase=True) for displ_cart in ndispl_carts
            ]
        else:
            ndispl_carts = displ_carts
            nevecs = [
                displacement_cart_to_evec(
                    displ_cart,
                    masses=masses,
                    scaled_positions=scaled_positions,
                    qpoint=qpt,
                    add_phase=True) for displ_cart in ndispl_carts
            ]

        for i in range(nbranch):
            freq = freqs[i]

            if freq < 0:
                evec = nevecs[i]
                #displ_cart=displ_carts[i]
                #evec = displacement_cart_to_evec(
                #    displ_cart,
                #    masses=masses,
                #    scaled_positions=scaled_positions,
                #    qpoint=qpt,
                #    add_phase=True)
                #evec=np.real(evec/np.linalg.norm(evec))

                scell = distorted_cell(atoms, supercell_matrix=np.eye(3) * 2)
                disp = scell._get_displacements(
                    evec, qpt, amplitude=1, argument=0)
                #disp=
                newcell = scell._get_cell_with_modulation(disp)
                newcell = Atoms(newcell)
                spacegroup = spglib.get_spacegroup(newcell, symprec=1e-3)
                print(spacegroup)

                write('unstable_modes/%s_%s.cif' % (qname, i), newcell)
                myfile.write("==============================\n")
                #evec=align(evec)
                #evec=np.real(evec/np.linalg.norm(evec))
                #evec=align([evec])[0]
                myfile.write(
                    "Evec: %s\n" % np.real(evec.reshape((16, 3))))
                myfile.write("Freq: %s\n" % (freq))
                myfile.write("%s_%s : %s\n" % (qname, i, spacegroup))

            #print(evec / np.linalg.norm(evec))
    return phbst


read_modes()
