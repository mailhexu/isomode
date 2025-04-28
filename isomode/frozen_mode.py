import numpy as np
import sys
from ase import Atoms
from ase.io import write
import numpy as np
import copy
import isomode.perovskite_mode as perovskite_mode
import spglib.spglib
from isomode.extern.cells import Supercell, Primitive

def default_mod_func(R):
    return 1.0

class distorted_cell():
    def __init__(self, atoms, supercell_matrix=np.eye(3)):
        self._primitive_cell = atoms
        self._supercell_matrix = supercell_matrix
        self._supercell = get_supercell(atoms, supercell_matrix, symprec=1e-4)
        self._N = np.linalg.det(supercell_matrix)

    def _get_cell_with_modulation(self, modulation):
        """
        x=x+dx
        """
        lattice = copy.deepcopy(self._supercell.get_cell())
        positions = copy.deepcopy(self._supercell.get_positions())
        masses = self._supercell.get_masses()
        symbols = self._supercell.get_chemical_symbols()
        positions += modulation.real
        scaled_positions = np.dot(positions, np.linalg.inv(lattice))
        for p in scaled_positions:
            p -= np.floor(p)
            cell = copy.deepcopy(self._supercell)
            cell.set_scaled_positions(scaled_positions)
        return cell

    def _get_displacements(self,
                           eigvec,
                           q,
                           amplitude,
                           argument,
                           mod_func=default_mod_func,
                           use_isotropy_amplitude=True):
        """
        displacements from eigvec, q, amplitude
        """
        m = self._supercell.get_masses()
        s2u_map = self._supercell.get_supercell_to_unitcell_map()
        u2u_map = self._supercell.get_unitcell_to_unitcell_map()
        s2uu_map = [u2u_map[x] for x in s2u_map]
        spos = self._supercell.get_scaled_positions()
        dim = self._supercell.get_supercell_matrix()
        r=np.dot(spos, dim.T)
        coefs = np.exp(2j * np.pi * np.dot(
            r, q)) * mod_func(r)  # Here Do not use sqrt(m)/ np.sqrt(m)
        u = []
        for i, coef in enumerate(coefs):
            eig_index = s2uu_map[i] * 3
            u.append(eigvec[eig_index:eig_index + 3] * coef)

        #u = np.array(u) / np.sqrt(len(m))
        u = np.array(u) / np.linalg.norm(u)  #/np.sqrt(self._N)
        phase_factor = self._get_phase_factor(u, argument)

        if use_isotropy_amplitude:
            amplitude = amplitude  #*self._N
        u *= phase_factor * amplitude

        return u

    def _get_phase_factor(self, modulation, argument):
        u = np.ravel(modulation)
        index_max_elem = np.argmax(abs(u))
        max_elem = u[index_max_elem]
        phase_for_zero = max_elem / abs(max_elem)
        phase_factor = np.exp(1j * np.pi * argument / 180) / phase_for_zero
        return phase_factor



def get_supercell(unitcell, supercell_matrix, symprec=1e-5):
    return Supercell(unitcell, supercell_matrix, symprec=symprec)


def get_primitive(supercell, primitive_frame, symprec=1e-5):
    return Primitive(supercell, primitive_frame, symprec=symprec)


def trim_cell(relative_axes, cell, symprec):
    """
    relative_axes: relative axes to supercell axes
    Trim positions outside relative axes

    """
    positions = cell.get_scaled_positions()
    numbers = cell.get_atomic_numbers()
    masses = cell.get_masses()
    #magmoms = cell.get_magnetic_moments()
    lattice = cell.get_cell()
    trimed_lattice = np.dot(relative_axes.T, lattice)

    trimed_positions = []
    trimed_numbers = []
    if masses is None:
        trimed_masses = None
    else:
        trimed_masses = []
    #if magmoms is None:
    #    trimed_magmoms = None
    #else:
    #    trimed_magmoms = []
    extracted_atoms = []

    positions_in_new_lattice = np.dot(positions,
                                      np.linalg.inv(relative_axes).T)
    positions_in_new_lattice -= np.floor(positions_in_new_lattice)
    trimed_positions = np.zeros_like(positions_in_new_lattice)
    num_atom = 0

    mapping_table = np.arange(len(positions), dtype='intc')
    symprec2 = symprec**2
    for i, pos in enumerate(positions_in_new_lattice):
        is_overlap = False
        if num_atom > 0:
            diff = trimed_positions[:num_atom] - pos
            diff -= np.rint(diff)
            # Older numpy doesn't support axis argument.
            # distances = np.linalg.norm(np.dot(diff, trimed_lattice), axis=1)
            # overlap_indices = np.where(distances < symprec)[0]
            distances2 = np.sum(np.dot(diff, trimed_lattice)**2, axis=1)
            overlap_indices = np.where(distances2 < symprec2)[0]
            if len(overlap_indices) > 0:
                is_overlap = True
                mapping_table[i] = extracted_atoms[overlap_indices[0]]

        if not is_overlap:
            trimed_positions[num_atom] = pos
            num_atom += 1
            trimed_numbers.append(numbers[i])
            if masses is not None:
                trimed_masses.append(masses[i])
            #if magmoms is not None:
            #    trimed_magmoms.append(magmoms[i])
            extracted_atoms.append(i)

    trimed_cell = Atoms(
        numbers=trimed_numbers,
        masses=trimed_masses,
        #magmoms=trimed_magmoms,
        scaled_positions=trimed_positions[:num_atom],
        cell=trimed_lattice, )

    return trimed_cell, extracted_atoms, mapping_table


def print_cell(cell, mapping=None, stars=None):
    symbols = cell.get_chemical_symbols()
    masses = cell.get_masses()
    magmoms = cell.get_magnetic_moments()
    lattice = cell.get_cell()
    print("Lattice vectors:")
    print("  a %20.15f %20.15f %20.15f" % tuple(lattice[0]))
    print("  b %20.15f %20.15f %20.15f" % tuple(lattice[1]))
    print("  c %20.15f %20.15f %20.15f" % tuple(lattice[2]))
    print("Atomic positions (fractional):")
    for i, v in enumerate(cell.get_scaled_positions()):
        num = " "
        if stars is not None:
            if i in stars:
                num = "*"
        num += "%d" % (i + 1)
        line = ("%5s %-2s%18.14f%18.14f%18.14f" %
                (num, symbols[i], v[0], v[1], v[2]))
        if masses is not None:
            line += " %7.3f" % masses[i]
        if magmoms is not None:
            line += "  %5.3f" % magmoms[i]
        if mapping is None:
            print(line)
        else:
            print(line + " > %d" % (mapping[i] + 1))



if __name__ == '__main__':
    test()

