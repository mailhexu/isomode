import pickle
from pyDFTutils.ase_utils import vesta_view, substitute_atoms
from ase.atoms import Atoms


def substitude_atoms(atoms, A1='Cs', A2='La', vca_format='abinit'):
    """
    for vasp Ti/Nb -> Ne/Nb
    for abinit Ti/Nb -> He
    """
    symbols = atoms.get_chemical_symbols()
    pos = atoms.get_scaled_positions()
    natoms = len(atoms)
    cell = atoms.get_cell()
    newsymbols = ''
    newpos = []
    for i in range(natoms):
        if symbols[i] == 'K':
            newsymbols+=A1
            newpos.append(pos[i])
        elif symbols[i] == 'Pb':
            newsymbols+=A2
            newpos.append(pos[i])
        elif symbols[i] == 'Nb':
            newsymbols+='Ti'
            newpos.append(pos[i])
        elif symbols[i] == 'Ta' and vca_format=='vasp':
            newsymbols+='He'
            newsymbols+='Nb'
            newpos.append(pos[i])
            newpos.append(pos[i])
        elif symbols[i] == 'Ta' and vca_format=='abinit':
            newsymbols+='He'
            newpos.append(pos[i])
    newatoms = Atoms(symbols=newsymbols, scaled_positions=newpos, cell=cell)
    return newatoms


def gen_all_structures(
        A1='Cs',
        A2='La', ):
    with open('single_mode.pickle', 'rb') as myfile:
        structures = pickle.load(myfile)
    names = []
    distorted_structures = []
    for structure in structures:
        label, direction, sym, _, atoms = structure
        name = "%s_%s" % (label, direction)
        names.append(name)
        newatoms = substitude_atoms(atoms, A1=A1, A2=A2)
        distorted_structures.append(newatoms)
    return names, distorted_structures

print(gen_all_structures())
