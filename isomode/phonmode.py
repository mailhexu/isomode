import os
import numpy as np
from scipy.optimize import minimize_scalar
from abipy.abilab import abiopen
from isomode.frozen_mode import distorted_cell
from ase import Atoms
from ase.io import write
import spglib.spglib
import matplotlib.pyplot as plt
from functools import partial
from numpy import array
import pickle


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
    if length == 3:
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
    theta = res.x
    new_modes = rotate_evecs(modes, theta)
    #print("theta: ", theta / np.pi)
    for i in range(len(modes)):
        mi = modes[i] / np.linalg.norm(modes[i])
        nmi = new_modes[i] / np.linalg.norm(new_modes[i])
        #print("Aligning=========\n", np.real(mi).reshape((16,3)), np.real(nmi).reshape((16,3)))
    return new_modes


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
                if len(multi_modes) == 3:  # multi, align and save
                    print(evals[i - 1], evec)
            # clean multi_modes and save new
            multi_modes = []
            multi_modes.append(evec)
        else:
            multi_modes.append(evec)
    return new_evecs


class phonon_distort_generator(object):
    def __init__(self, fname, qdict):
        self.fname = fname
        self.qdict = qdict
        self.read_phbst_file()
        self.distorted_structures = []

    def read_phbst_file(self):
        self.phbst = abiopen(self.fname)
        self.atoms = self.phbst.structure.to_ase_atoms()
        # write('primitive.cif',atoms)
        self.nbranch = 3 * len(self.atoms)
        self.masses = self.atoms.get_masses()
        self.scaled_positions = self.atoms.get_scaled_positions()

    def generate_primitive_structure(self, path=None):
        self.atoms.set_pbc([1, 1, 1])
        if path is not None:
            if not os.path.exists(path):
                os.makedirs(path)
            write(os.path.join(path, 'primitive.cif'), self.atoms)
        return self.atoms

    def generate_combined_mode_distorted_structures(
            self,
            supercell_matrix=np.eye(3) * 2,
            amplitude=1.0,
            max_freq=0.0,
            path=None,
            align_evecs=False,
            align_disp=False, ):
        """
        qname, qpoint, ibranch, freq, structure, spacegroup
        """
        if not os.path.exists(path):
            os.makedirs(path)

        disp=0.0
        for qname, qpt in self.qdict.items():
            displ_carts = [
                self.phbst.get_phmode(qpt, i).displ_cart
                for i in range(self.nbranch)
            ]
            freqs = [
                self.phbst.get_phmode(qpt, i).freq
                for i in range(self.nbranch)
            ]
            if align_evecs:
                evecs = [
                    displacement_cart_to_evec(
                        displ_cart,
                        masses=self.masses,
                        scaled_positions=self.scaled_positions,
                        qpoint=qpt,
                        add_phase=True) for displ_cart in displ_carts
                ]

                nevecs = align_all_modes(freqs, evecs)
            elif align_disp:
                ndispl_carts = align_all_modes(freqs, displ_carts)
                nevecs = [
                    displacement_cart_to_evec(
                        displ_cart,
                        masses=self.masses,
                        scaled_positions=self.scaled_positions,
                        qpoint=qpt,
                        add_phase=True) for displ_cart in ndispl_carts
                ]
            else:
                ndispl_carts = displ_carts
                nevecs = [
                    displacement_cart_to_evec(
                        displ_cart,
                        masses=self.masses,
                        scaled_positions=self.scaled_positions,
                        qpoint=qpt,
                        add_phase=True) for displ_cart in ndispl_carts
                ]
            scell = distorted_cell(
                        self.atoms, supercell_matrix=supercell_matrix)
            for i in range(self.nbranch):
                freq = freqs[i]
                if freq < max_freq:
                    evec = nevecs[i]
                    scell = distorted_cell(
                        self.atoms, supercell_matrix=supercell_matrix)
                    disp += scell._get_displacements(
                        evec, qpt, amplitude=amplitude, argument=0)
            newcell = scell._get_cell_with_modulation(disp)
            newcell = Atoms(newcell)
            spacegroup = spglib.get_spacegroup(newcell, symprec=1e-3)
            if path is not None:
                newcell.set_pbc([1, 1, 1])
                write(os.path.join(path, 'combined_modes.cif'), newcell)

    def generate_distorted_structures(self,
                                      supercell_matrix=np.eye(3) * 2,
                                      amplitude=1.0,
                                      unstable_only=True,
                                      path=None):
        """
        qname, qpoint, ibranch, freq, structure, spacegroup
        """
        if not os.path.exists(path):
            os.makedirs(path)
        for qname, qpt in self.qdict.items():
            displ_carts = [
                self.phbst.get_phmode(qpt, i).displ_cart
                for i in range(self.nbranch)
            ]
            freqs = [
                self.phbst.get_phmode(qpt, i).freq
                for i in range(self.nbranch)
            ]
            align_evecs = False
            align_disp = False
            if align_evecs:
                evecs = [
                    displacement_cart_to_evec(
                        displ_cart,
                        masses=self.masses,
                        scaled_positions=self.scaled_positions,
                        qpoint=qpt,
                        add_phase=True) for displ_cart in displ_carts
                ]

                nevecs = align_all_modes(freqs, evecs)
            elif align_disp:
                ndispl_carts = align_all_modes(freqs, displ_carts)
                nevecs = [
                    displacement_cart_to_evec(
                        displ_cart,
                        masses=self.masses,
                        scaled_positions=self.scaled_positions,
                        qpoint=qpt,
                        add_phase=True) for displ_cart in ndispl_carts
                ]
            else:
                ndispl_carts = displ_carts
                nevecs = [
                    displacement_cart_to_evec(
                        displ_cart,
                        masses=self.masses,
                        scaled_positions=self.scaled_positions,
                        qpoint=qpt,
                        add_phase=True) for displ_cart in ndispl_carts
                ]

            for i in range(self.nbranch):
                freq = freqs[i]
                if freq < 0 or (not unstable_only):
                    print(len(nevecs))
                    evec = displ_carts[i]#nevecs[i]
                    print(evec)
                    scell = distorted_cell(
                        self.atoms, supercell_matrix=supercell_matrix)
                    disp = scell._get_displacements(
                        evec, qpt, amplitude=amplitude, argument=0)
                    newcell = scell._get_cell_with_modulation(disp)
                    newcell = Atoms(newcell)
                    spacegroup = spglib.get_spacegroup(newcell, symprec=1e-3)
                    if path is not None:
                        newcell.set_pbc([1, 1, 1])
                        write(
                            os.path.join(path, '%s_%s.cif' % (qname, i)),
                            newcell)
                    self.distorted_structures.append({
                        'qname': qname,
                        'qpoint': qpt,
                        'ibranch': i,
                        'evec': evec,
                        'amplitude': amplitude,
                        'spacegroup': spacegroup,
                        'atoms': newcell
                    })
        if path is not None:
            pfname = os.path.join(path, 'phmodes.pickle')
            with open(pfname, 'wb') as myfile:
                pickle.dump(self.distorted_structures, myfile)
        return self.distorted_structures


def test_phbst_modes(fname='../data/DionJ/phonon/run.abo_PHBST.nc'):
    qdict = {
        'Gamma': [0.0, 0.0, 0.0],
        'Xy': [0, 0.5, 0],
        #'Xx': [0.5, 0.0, 0],
        'M': [0.5, 0.5, 0],
        'Rx': [0.5, 0.0, 0.5],
        #'Ry': [0.0, 0.5, 0.5],
        'A': [0.5, 0.5, 0.5],
        'Z': [0, 0, 0.5]
    }
    dg = phonon_distort_generator(fname, qdict)
    dg.generate_primitive_structure(path='nmodes')
    ds = dg.generate_distorted_structures(
        supercell_matrix=np.eye(3) * 2,
        amplitude=0.4,
        unstable_only=True,
        path='nmodes')
    print(ds)


#test_phbst_modes()
