import numpy as np
import sys
from ase import Atoms
from ase.io import write
import numpy as np
import copy
import isomode.perovskite_mode as perovskite_mode
import spglib.spglib

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


class Supercell(Atoms):
    """Build supercell from supercell matrix
    In this function, unit cell is considered
    [1,0,0]
    [0,1,0]
    [0,0,1].
    Supercell matrix is given by relative ratio, e.g,
    [-1, 1, 1]
    [ 1,-1, 1]  is for FCC from simple cubic.
    [ 1, 1,-1]
    In this case multiplicities of surrounding simple lattice are [2,2,2].

    First, create supercell with surrounding simple lattice.
    Second, trim the surrounding supercell with the target lattice.
    """

    def __init__(self, unitcell, supercell_matrix, symprec=1e-5):
        self._s2u_map = None
        self._u2s_map = None
        self._u2u_map = None
        self._supercell_matrix = np.array(supercell_matrix, dtype='intc')
        self._create_supercell(unitcell, symprec)

    def get_supercell_matrix(self):
        return self._supercell_matrix

    def get_supercell_to_unitcell_map(self):
        return self._s2u_map

    def get_unitcell_to_supercell_map(self):
        return self._u2s_map

    def get_unitcell_to_unitcell_map(self):
        return self._u2u_map

    def _create_supercell(self, unitcell, symprec):
        mat = self._supercell_matrix
        frame = self._get_surrounding_frame(mat)
        sur_cell, u2sur_map = self._get_simple_supercell(frame, unitcell)

        # Trim the simple supercell by the supercell matrix
        trim_frame = np.array([
            mat[0] / float(frame[0]), mat[1] / float(frame[1]),
            mat[2] / float(frame[2])
        ])
        supercell, sur2s_map, mapping_table = trim_cell(trim_frame, sur_cell,
                                                        symprec)

        multi = supercell.get_number_of_atoms(
        ) // unitcell.get_number_of_atoms()

        if multi != determinant(self._supercell_matrix):
            print("Supercell creation failed.")
            print("Probably some atoms are overwrapped. "
                  "The mapping table is give below.")
            print(mapping_table)
            Atoms.__init__(self)
        else:
            Atoms.__init__(
                self,
                numbers=supercell.get_atomic_numbers(),
                masses=supercell.get_masses(),
                #magmoms=supercell.get_magnetic_moments(),
                scaled_positions=supercell.get_scaled_positions(),
                cell=supercell.get_cell(),
                #pbc=True)
            )
            self._u2s_map = np.arange(unitcell.get_number_of_atoms()) * multi
            self._u2u_map = dict([(j, i) for i, j in enumerate(self._u2s_map)])
            self._s2u_map = np.array(u2sur_map)[sur2s_map] * multi

    def _get_surrounding_frame(self, supercell_matrix):
        # Build a frame surrounding supercell lattice
        # For example,
        #  [2,0,0]
        #  [0,2,0] is the frame for FCC from simple cubic.
        #  [0,0,2]
        m = np.array(supercell_matrix)
        axes = np.array([[0, 0, 0], m[:, 0], m[:, 1], m[:, 2],
                         m[:, 1] + m[:, 2], m[:, 2] + m[:, 0],
                         m[:, 0] + m[:, 1], m[:, 0] + m[:, 1] + m[:, 2]])
        frame = [max(axes[:, i]) - min(axes[:, i]) for i in (0, 1, 2)]
        return frame

    def _get_simple_supercell(self, multi, unitcell):
        # Scaled positions within the frame, i.e., create a supercell that
        # is made simply to multiply the input cell.
        positions = unitcell.get_scaled_positions()
        numbers = unitcell.get_atomic_numbers()
        masses = unitcell.get_masses()
        #magmoms = unitcell.get_magnetic_moments()
        lattice = unitcell.get_cell()

        atom_map = []
        positions_multi = []
        numbers_multi = []
        if masses is None:
            masses_multi = None
        else:
            masses_multi = []
        #if magmoms is None:
        #    magmoms_multi = None
        #else:
        #    magmoms_multi = []
        for l, pos in enumerate(positions):
            for i in range(multi[2]):
                for j in range(multi[1]):
                    for k in range(multi[0]):
                        positions_multi.append([(pos[0] + k) / multi[0],
                                                (pos[1] + j) / multi[1],
                                                (pos[2] + i) / multi[2]])
                        numbers_multi.append(numbers[l])
                        if masses is not None:
                            masses_multi.append(masses[l])
                        atom_map.append(l)
                        #if magmoms is not None:
                        #    magmoms_multi.append(magmoms[l])

        simple_supercell = Atoms(
            numbers=numbers_multi,
            masses=masses_multi,
            #magmoms=magmoms_multi,
            scaled_positions=positions_multi,
            cell=np.dot(np.diag(multi), lattice),
            #pbc=True
        )

        return simple_supercell, atom_map


class Primitive(Atoms):
    def __init__(self, supercell, primitive_matrix, symprec=1e-5):
        """
        primitive_matrix (3x3 matrix):
        Primitive lattice is given with respect to supercell by
           np.dot(primitive_matrix.T, supercell.get_cell())
        """
        self._primitive_matrix = np.array(primitive_matrix)
        self._symprec = symprec
        self._p2s_map = None
        self._s2p_map = None
        self._p2p_map = None
        self._primitive_cell(supercell)
        self._supercell_to_primitive_map(supercell.get_scaled_positions())
        self._primitive_to_primitive_map()

    def get_primitive_matrix(self):
        return self._primitive_matrix

    def get_primitive_to_supercell_map(self):
        return self._p2s_map

    def get_supercell_to_primitive_map(self):
        return self._s2p_map

    def get_primitive_to_primitive_map(self):
        return self._p2p_map

    def _primitive_cell(self, supercell):
        trimed_cell, p2s_map, mapping_table = trim_cell(
            self._primitive_matrix, supercell, self._symprec)
        Atoms.__init__(
            self,
            numbers=trimed_cell.get_atomic_numbers(),
            masses=trimed_cell.get_masses(),
            #magmoms=trimed_cell.get_magnetic_moments(),
            scaled_positions=trimed_cell.get_scaled_positions(),
            cell=trimed_cell.get_cell(),
            pbc=True)

        self._p2s_map = np.array(p2s_map, dtype='intc')

    def _supercell_to_primitive_map(self, pos):
        inv_F = np.linalg.inv(self._primitive_matrix)
        s2p_map = []
        for i in range(pos.shape[0]):
            s_pos = np.dot(pos[i], inv_F.T)
            for j in self._p2s_map:
                p_pos = np.dot(pos[j], inv_F.T)
                diff = p_pos - s_pos
                diff -= np.rint(diff)
                if (abs(diff) < self._symprec).all():
                    s2p_map.append(j)
                    break
        self._s2p_map = np.array(s2p_map, dtype='intc')

    def _primitive_to_primitive_map(self):
        """
        Mapping table from supercell index to primitive index
        in primitive cell
        """
        self._p2p_map = dict([(j, i) for i, j in enumerate(self._p2s_map)])


def determinant(m):
    return (m[0][0] * m[1][1] * m[2][2] - m[0][0] * m[1][2] * m[2][1] + m[0][1]
            * m[1][2] * m[2][0] - m[0][1] * m[1][0] * m[2][2] + m[0][2] *
            m[1][0] * m[2][1] - m[0][2] * m[1][1] * m[2][0])

if __name__ == '__main__':
    test()

