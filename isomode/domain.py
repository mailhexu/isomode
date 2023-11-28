#/usr/bin/env python
"""
Generate domains and domain walls from phonon.
"""

from .frozen_mode import distorted_cell, default_mod_func
from pyDFTutils.ase_utils import vesta_view
from dataclasses import dataclass



class DomainBuilder():
    def __init__(self, primcell, supercell_matrix, mode_list=[]):
        """
        primitive cell: primitive cell.
        displacement mode list: [displacement vector, amplitude, argument, modulation_function]
        """
        self.primcell = primcell
        self.mode_list = mode_list
        self.supercell_matrix = supercell_matrix
        self.supercell = distorted_cell(
            self.primcell, supercell_matrix=self.supercell_matrix)

    def add_mode(self,
                 displacement,
                 qpoint,
                 amplitude,
                 argument=0,
                 modulation_function=default_mod_func):
        self.mode_list.append(
            [displacement, qpoint, amplitude, argument, modulation_function])

    def get_distortions(self):
        distortion=0.0
        for disp, q, amp, arg, mod_func in self.mode_list:
            distortion = self.supercell._get_displacements(
                self,
                disp,
                q,
                amp,
                arg,
                mod_func=mod_func,
                use_isotropy_amplitude=True)
            distortion += np.array(distortion)
        return distortion

    def get_supercell(self):
        return self.supercell._get_cell_with_modulation(self.get_distortions())

    def view_supercell(self):
        vesta_view(self.get_supercell())

    def save_supercell(self):
        pass



class DomainBuilderPhonon(DomainBuilder):
    def __init__(self, phonon, supercell_matrix, mode_list=[]):
        self.supercell=phonon.supercell

    def add_mode(self, phonon, imode, iq, amplitude, argument=0, modulation_function):
        pass
