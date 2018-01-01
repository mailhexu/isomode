import numpy as np
import os
import spglib
from ase import Atoms
from ase.io import read, write
from pyDFTutils.ase_utils import symbol_number, symnum_to_sym, vesta_view
from collections import OrderedDict, defaultdict
import re
from itertools import combinations
import copy
import json


def parse_mode_name_line(line=None):
    """
    return a dictionary ret
    ret['kpt']: the kpoint. a list. Note: In python2, 1/2 -> 0 . So it could be wrong in python2.
    ret['kpt_string']: kpt in string. 1/2 is still 1/2
    ret['normfactor']
    ret['label']: eg GM3+
    ret['symmetry']: eg A2u
    fullname:  '[0,0,0]GM1+[Nb2:g:dsp]A1(a)'
    direction: a/b/c
    """
    kpt_string = re.findall(r'\[.*\d\]', line)[0]
    kpt = (eval(kpt_string))
    normfactor_string = re.findall(r'normfactor\s*=\s*(.*)', line)[0]
    normfactor = float(normfactor_string)
    label, symmetry = re.findall(r'\](.*?)\(', line)

    a1 = re.findall(r'(\[.*?)\(', line)[0]
    a2 = re.findall(r'\[.*?\)', line)[1]
    fullname = a1 + a2
    direction = fullname[-2]

    return {
        'kpt': kpt,
        'kpt_string': kpt_string,
        'normfactor': normfactor,
        'label': label,
        'symmetry': symmetry,
        'direction': direction,
        'fullname': fullname,
    }


def parse_mode_amp_line(line):
    """
    return dict, of which the keys are
    kpt
    kpt_string
    label
    symmetry
    direction
    fullname
    amp
    """
    fullname, amp, _, _ = line.strip().split()
    amp = float(amp)
    kpt_string = re.findall(r'\[.*\d\]', line)[0]
    kpt = (eval(kpt_string))
    label = re.findall(r'\](.*?)\[', line)[0]
    symmetry = re.findall(r'.*\](.*?)\(',
                          line)[0]  # first .* is greedy to eat all before ]
    direction = fullname[-2]

    return {
        'kpt': kpt,
        'kpt_string': kpt_string,
        'label': label,
        'symmetry': symmetry,
        'direction': direction,
        'fullname': fullname,
        'amp': amp
    }


#print(parse_mode_amp_line(
#    "[0,0,0]GM5-[Cs1:d:dsp]Eu(a)         0.00000   0.00000  0.00000"))


def find_primitive(atoms, symprec=1e-4):
    """
    find the primitive cell withh regard to the magnetic structure. a atoms object is returned.
    """
    #atoms_mag,sym_dict=ref_atoms_mag(atoms)
    cell, scaled_pos, chem_nums = spglib.find_primitive(atoms, symprec=symprec)
    #print(spglib.get_symmetry(atoms, symprec=1e-5))
    chem_sym = 'H%d' % (len(chem_nums))
    new_atoms = Atoms(chem_sym)

    new_atoms.set_atomic_numbers(chem_nums)
    new_atoms.set_cell(cell)
    new_atoms.set_scaled_positions(scaled_pos)
    #new_atoms=rev_ref_atoms(new_atoms,sym_dict)
    return new_atoms


#print(
#    parse_mode_name_line(
#        line="P4/mmm[0,0,0]GM1+(a)[Nb2:g:dsp]A1(a) normfactor = 0.00825"))


def split_symnum(symnum):
    """
    symnum-> sym. eg: Fe1-> Fe
    """
    try:
        a = re.search('[A-Za-z]+', symnum).group()
        b = int(symnum[len(a):])
        return a, b
    except AttributeError:
        raise AttributeError('%s is not a good symbol_number' % symnum)


def read_mode_amplitudes(filename='mtable.org',
                         use_direction=True,
                         randomize=False,
                         use_abs=False,
                         use_max=False):
    detail_mode_amplitude_dict = defaultdict(float)
    total_mode_dict = defaultdict(float)
    total_mode_definitions = dict()
    with open(filename) as myfile:
        lines = myfile.readlines()
    for line in lines:
        try:
            result = parse_mode_amp_line(line)
        except Exception:
            result = None
        if result is not None:
            fullname = result['fullname']
            label = result['label']
            direction = result['direction']
            if use_direction:
                shortname = (label, direction)
            else:
                shortname = label

            amp = result['amp']
            if use_abs:
                amp = abs(amp)
            if np.abs(amp) > 0.001:
                if randomize:
                    amp = amp * (1.0 + (random.random() - 1.0) * 0.03)
                detail_mode_amplitude_dict[fullname] += amp
                total_mode_dict[shortname] += amp
                if shortname not in total_mode_definitions:
                    total_mode_definitions[shortname] = defaultdict(float)
                else:
                    total_mode_definitions[shortname][fullname] += amp
            #print(m['label'], m['direction'], m['amp'])
    if use_max:
        total_mode_definitions_max = dict()
        for key, val in total_mode_definitions.items():
            #print(max(val.items(), key=lambda x: x[1]))
            total_mode_definitions_max[key] = dict(
                [max(val.items(), key=lambda x: abs(x[1]))])
        total_mode_definitions = total_mode_definitions_max
    return detail_mode_amplitude_dict, total_mode_dict, total_mode_definitions


def test_read_mode_amplitudes():
    print(read_mode_amplitudes('mtable.org', use_abs=True, use_max=True)[2])


#print(read_mode_amplitudes('mtable.org', use_direction=False, use_abs=True, use_max=False)[2])


class Isomode(object):
    def __init__(self, fname):
        self.fname = fname
        with open(self.fname) as myfile:
            self.lines = myfile.readlines()
        self.atoms = None
        self.cellpars = []
        self.symbols = []
        self.positions = []
        self.symdict = {}
        self.natom = None
        self.mode_definitions = {}
        self.read_supercell()
        self.read_mode_definitions()

    def read_primitive_cell(self):
        pass

    def read_supercell(self):
        inside = False
        sympos = []
        iatom = 0
        for iline, line in enumerate(self.lines):
            if line.strip().startswith('Undistorted superstructure'):
                inside = True
                continue
            if inside:
                if line.strip().startswith('a='):
                    # eg. a=7.80962, b=7.80962, c=30.28940, alpha=90.00000, beta=90.00000, gamma=90.00000
                    segs = line.strip().split(',')
                    for seg in segs:
                        self.cellpars.append(float(seg.strip().split('=')[1]))
                elif line.strip().startswith('atom'):
                    pass
                elif line.strip() == '':
                    inside = False
                else:
                    symnum, site, x, y, z, occ, displ = line.strip().split()
                    sym, num = split_symnum(symnum)
                    x = float(x)
                    y = float(y)
                    z = float(z)
                    occ = float(occ)
                    displ = float(displ)
                    sympos.append([symnum, x, y, z])
                    self.symbols.append(sym)
                    self.positions.append(np.array([x, y, z]))
                    self.symdict[symnum] = iatom
                    iatom = iatom + 1

        self.atoms = Atoms(
            symbols=self.symbols,
            scaled_positions=self.positions,
            cell=self.cellpars)
        self.natom = len(self.atoms)

    def substitute_element(self, subs_dict):
        symbols=self.symbols
        new_symbols=[]
        for sym in symbols:
            if sym in subs_dict:
                new_symbols.append(subs_dict[sym])
            else:
                new_symbols.append(sym)
        self.symbols=new_symbols
        self.atoms = Atoms(
            symbols=self.symbols,
            scaled_positions=self.positions,
            cell=self.cellpars)


    def read_mode_definitions(self):
        """
        read displacive mode definitions
        """
        inside = False
        mode_definitions = OrderedDict()
        inside_mode = False
        for iline, line in enumerate(self.lines):
            if line.strip().startswith("Displacive mode definitions"):
                inside = True
                continue
            elif line.strip().startswith("Displacive mode amplitudes"):
                inside = False
            elif inside:
                if line.find('normfactor') != -1:  # begining of one mode
                    nameline = line
                    r = parse_mode_name_line(line=nameline)
                    mode_info = {}
                    for key in r:
                        mode_info[key] = r[key]
                    kpt = r['kpt']
                    kpt_string = r['kpt_string']
                    normfactor = r['normfactor']
                    label = r['label']
                    symmetry = r['symmetry']
                    fullname = r['fullname']
                    inside_mode = True
                    #deltas = {}
                    mode = np.zeros([self.natom, 3], dtype=float)
                    continue
                if inside_mode:
                    if line.strip() == '':  # end of one mode.
                        inside_mode = False
                        self.mode_definitions[mode_info['fullname']] = {
                            'mode_info': mode_info,
                            'mode': mode
                        }
                    elif line.strip().startswith('atom'):
                        pass
                    else:  # a line  of displacement.
                        symnum, x, y, z, dx, dy, dz = line.strip().split()
                        #sym, num = split_symnum(symnum)
                        _, _, _, dx, dy, dz = map(float, (x, y, z, dx, dy, dz))
                        delta = (dx, dy, dz)
                        #deltas[symnum] = delta
                        mode[self.symdict[symnum]] = delta

    def get_mode_info(self, fullname, info_key):
        """
        return information of a mode by fullname
        fullname: eg.  [0,0,0]GM5-[Cs1:d:dsp]Eu(a)
        info_key: kpt|kpt_string|norm_factor|label|symmetry|fullanme|delta
        """
        if (info_key == 'displacement') or (info_key == 'mode'):
            return self.mode_definitions[fullname]['mode']
        else:
            return self.mode_definitions[fullname]['mode_info'][info_key]

    def get_mode_displacement(self, fullname):
        return self.mode_definitions[fullname]['mode']

    def get_all_mode_fullname(self):
        return self.mode_definitions.keys()

    def get_undistorted_structure(self):
        return self.atoms

    def get_total_mode_displacement(self):
        pass

    def get_total_mode_info(self):
        pass

    def get_distorted_structure(self, mode_dict, primitive=True):
        distorted_structure = self.atoms.copy()
        scaled_positions = distorted_structure.get_scaled_positions()
        for fullname in mode_dict:
            amp = mode_dict[fullname]
            normfactor = self.mode_definitions[fullname]['mode_info'][
                'normfactor']
            disp = [
                normfactor * amp * np.dot(self.atoms.get_cell(), m)
                for m in self.mode_definitions[fullname]['mode']
            ]
            scaled_positions += disp

        distorted_structure.set_scaled_positions(scaled_positions)
        #print(spglib.get_spacegroup(distorted_structure, symprec=1e-4))
        if primitive:
            distorted_structure = find_primitive(
                distorted_structure, symprec=1e-4)
            #print(distorted_structure)
        return distorted_structure

    def read_mode_amplitudes(self, fname=None, use_direction=True, **kwargs):
        """
        it read amplitudes from a file.
        saved information:
        1. each detailed mode amplitude
        2. each total mode (defined by label+direction) amplitude
        3. a map to total_mode from detailed mode.
        """
        if fname is None:
            fname = self.fname
        (self.detail_mode_amplitude_dict, self.total_mode_dict,
         self.total_mode_definitions) = read_mode_amplitudes(
             fname, use_max=True)
        #print(self.total_mode_definitions)
        #print(self.total_mode_dict)

    def prepare_structure_single(self,
                                 json_file_name=None,
                                 amp=None,
                                 cif_dir=None,
                                 primitive=False):
        #print(self.total_mode_definitions)
        structure_list = []
        if cif_dir is not None:
            if not os.path.exists(cif_dir):
                os.makedirs(cif_dir)
        for key, mode_dict in self.total_mode_definitions.items():
            label, direction = key
            # Note sometimes 'b' without 'a'
            if (direction == 'a') or (
                (label, 'a') not in self.total_mode_definitions):
                if amp is not None:
                    for fullname in mode_dict:
                        mode_dict[fullname] = amp
                distorted_structure = self.get_distorted_structure(
                    mode_dict, primitive=primitive)
                #print(distorted_structure)
                spgroup = spglib.get_spacegroup(
                    distorted_structure, symprec=1e-4)
                filename = '%s_%s.cif' % (label, direction)
                print('%s %s %s %s\n' % (label, direction, spgroup, filename))
                structure_list.append(
                    [label, direction, spgroup, filename, distorted_structure])
                if cif_dir is not None:
                    cif_fname = os.path.join(cif_dir, filename)
                    write(cif_fname, distorted_structure)
        if json_file_name is not None:
            json_file = open(json_file_name, 'w')
            json_list = tuple(item[:-1] for item in structure_list)
            json.dump(json_list, json_file, indent=4)
            json_file.close()
        return structure_list

def read_unstable_modes(fname='./mtable.org'):
    with open(fname) as myfile:
        lines = myfile.readlines()
    for line in lines:
        if line.strip().startswith('[') and not line.find('all') != -1:
            m = parse_mode_amp_line(line)


modes = {
    '[0,0,1/2]Z5-[Ca1:h:dsp]E(a)': 1.0,
    '[0,0,1/2]Z5-[Ca1:h:dsp]E(b)': 1.0,
    '[0,0,0]GM5+[Ca1:h:dsp]E(a)': 1.0,
    '[0,0,0]GM5+[Ca1:h:dsp]E(b)': 1.0,
    '[0,0,0]GM5-[Cs1:d:dsp]Eu(a)': 1.0,
    '[0,0,0]GM5-[Cs1:d:dsp]Eu(b)': 1.0,
    '[1/2,1/2,1/2]A1-[O5:i:dsp]B1(a)': 1.0,
    '[1/2,1/2,1/2]A2+[O1:f:dsp]B2u(a)': 1.0,
    '[1/2,1/2,1/2]A3+[O1:f:dsp]B3u(a)': 1.0,
    '[1/2,1/2,0]M3+[O1:f:dsp]B3u(a)': 1.0,
    '[1/2,1/2,0]M5+[O1:f:dsp]B1u(a)': 1.0,
    '[0,1/2,0]X3+[O3:g:dsp]E(a)': 1.0,
    '[0,1/2,0]X3+[O3:g:dsp]E(b)': 1.0,
    '[0,1/2,1/2]R3+[O3:g:dsp]E(a)': 1.0,
    '[0,1/2,1/2]R3+[O3:g:dsp]E(b)': 1.0,
}


def label_direction_modes(modes):
    ldmodes = {}
    for mode in modes:
        fullname = mode
        label = re.findall(r'\](.*?)\[', fullname)[0]
        symmetry = re.findall(
            r'.*\](.*?)\(',
            fullname)[0]  # first .* is greedy to eat all before ]
        direction = fullname[-2]
        ldmode = (label, direction)
        #ldmodes.append(ldmode)
        ldmodes[ldmode] = fullname

    return ldmodes


def group_direction(modes, direction='a'):
    ldmodes = label_direction_modes(modes)
    return dict([(x, ldmodes[x]) for x in ldmodes if x[1] == direction])


class Iso_unstable_modes(object):
    def __init__(self, all_isomodes, mode_list):
        self.all_isomodes = all_isomodes
        self.mode_list = mode_list
        self.label_direction_modes = label_direction_modes(modes)

    def write_table(self, fname):
        pass

    def get_high_symmetry_structure(self):
        pass

    def get_high_symmetry_structure_suprecells(self):
        pass

    def gen_structure_a(self, amplitude=0.02):
        ldmodes = self.label_direction_modes
        modes_a = dict([(x, ldmodes[x]) for x in ldmodes if x[1] == 'a'])
        for (label, direction), fullname in modes_a.items():
            distorted_structure = self.all_isomodes.get_distorted_structure(
                modes)
            spgroup = spglib.get_spacegroup(distorted_structure, symprec=1e-4)

    def write_structure_a(amplitude=0.02):
        #myparser = isodistort_parser('./isodistort_modes_a.txt')
        myfile = open('single_mode.org', 'w')
        myfile.write('|label|ndirection|spacegroup|\n')
        for name, modes in mode_groups.items():
            label = myparser.modes[list(modes.keys())[0]]['mode_info']['label']
            n_direction = len(modes)
            distorted_structure = myparser.get_distorted_structure(modes)

            spgroup = spglib.get_spacegroup(distorted_structure, symprec=1e-4)
            myfile.write('|%s|%s|%s|%s|\n' % (label, n_direction, spgroup,
                                              name))
        myfile.close()


def gen_structure_a(modes, amplitude=0.01):
    myparser = isodistort_parser('./isodistort_modes_a.txt')
    myfile = open('single_mode.org', 'w')
    myfile.write('|label|ndirection|spacegroup|\n')
    for name, modes in mode_groups.items():
        label = myparser.modes[list(modes.keys())[0]]['mode_info']['label']
        n_direction = len(modes)
        distorted_structure = myparser.get_distorted_structure(modes)

        spgroup = spglib.get_spacegroup(distorted_structure, symprec=1e-4)
        myfile.write('|%s|%s|%s|%s|\n' % (label, n_direction, spgroup, name))
    myfile.close()


def gen_structure_ab():
    pass


def gen_structure_aba():
    pass


def gen_structure_abab():
    pass


def test():
    myparser = isodistort_parser('./isodistort_modes.txt')
    #print(myparser.modes)
    #vesta_view(myparser.get_distorted_structure({'[0,0,1/2]Z5-[O9:g:dsp]E(a)': 1.0}))
    #vesta_view(

    mode_groups = dict(
        Z5_m_1={'[0,0,1/2]Z5-[Ca1:h:dsp]E(a)': 1.0, },
        Z5_m_2={
            '[0,0,1/2]Z5-[Ca1:h:dsp]E(a)': 1.0,
            '[0,0,1/2]Z5-[Ca1:h:dsp]E(b)': 1.0,
        },
        GM5_p_1={'[0,0,0]GM5+[Ca1:h:dsp]E(a)': 1.0, },
        GM5_p_2={
            '[0,0,0]GM5+[Ca1:h:dsp]E(a)': 1.0,
            '[0,0,0]GM5+[Ca1:h:dsp]E(b)': 1.0
        },
        GM5_m_1={'[0,0,0]GM5-[Cs1:d:dsp]Eu(a)': 1.0, },
        GM5_m_2={
            '[0,0,0]GM5-[Cs1:d:dsp]Eu(a)': 1.0,
            '[0,0,0]GM5-[Cs1:d:dsp]Eu(b)': 1.0,
        },
        A1_m_1={'[1/2,1/2,1/2]A1-[O5:i:dsp]B1(a)': 1.0},
        A2_p_1={'[1/2,1/2,1/2]A2+[O1:f:dsp]B2u(a)': 1.0, },
        A3_p_1={'[1/2,1/2,1/2]A3+[O1:f:dsp]B3u(a)': 1.0, },
        M3_p_1={'[1/2,1/2,0]M3+[O1:f:dsp]B3u(a)': 1.0, },
        M5_p_1={'[1/2,1/2,0]M5+[O1:f:dsp]B1u(a)': 1.0, },
        X3_p_1={'[0,1/2,0]X3+[O3:g:dsp]E(a)': 1.0, },
        X3_p_2={
            '[0,1/2,0]X3+[O3:g:dsp]E(a)': 1.0,
            '[0,1/2,0]X3+[O3:g:dsp]E(b)': 1.0,
        },
        R3_p_1={'[0,1/2,1/2]R3+[O3:g:dsp]E(a)': 1.0, },
        R3_p_2={
            '[0,1/2,1/2]R3+[O3:g:dsp]E(a)': 1.0,
            '[0,1/2,1/2]R3+[O3:g:dsp]E(b)': 1.0,
        })
    myfile = open('single_mode.org', 'w')
    myfile.write('|label|ndirection|spacegroup|\n')
    for name, modes in mode_groups.items():
        label = myparser.modes[list(modes.keys())[0]]['mode_info']['label']
        n_direction = len(modes)
        distorted_structure = myparser.get_distorted_structure(modes)

        spgroup = spglib.get_spacegroup(distorted_structure, symprec=1e-4)
        myfile.write('|%s|%s|%s|%s|\n' % (label, n_direction, spgroup, name))
    myfile.close()

    myfile = open('double_mode.org', 'w')
    myfile.write('|label1|ndirection1|label2|ndirection2|spacegroup|\n')
    for comb in combinations(mode_groups.items(), 2):
        g1, g2 = comb
        #print(g1)
        key, modes1 = g1
        label1 = myparser.modes[list(modes1.keys())[0]]['mode_info']['label']
        n_direction1 = len(modes1)

        key, modes2 = g2
        label2 = myparser.modes[list(modes2.keys())[0]]['mode_info']['label']
        n_direction2 = len(modes2)

        if label1 == label2:
            continue
        else:
            modes = copy.deepcopy(modes1)
            modes.update(modes2)
            distorted_structure = myparser.get_distorted_structure(modes)
            spgroup = spglib.get_spacegroup(distorted_structure, symprec=1e-4)
            myfile.write('|%s|%s|%s|%s|%s|\n' % (label1, n_direction1, label2,
                                                 n_direction2, spgroup))
    myfile.close()

    myparser.get_distorted_structure({
        '[0,0,1/2]Z5-[Ca1:h:dsp]E(a)': 1.0,
        #'[0,0,1/2]Z5-[Ca1:h:dsp]E(b)':1.0,
        '[0,0,1/2]Z5-[O1:f:dsp]B3u(a)': 1.0,
        #'[0,0,1/2]Z5-[O1:f:dsp]B3u(b)':1.0,
        '[0,0,1/2]Z5-[O9:g:dsp]E(a)': 1.0,
        #'[0,0,1/2]Z5-[Nb1:a:dsp]Eu(b)': 1.0,
        #'[0,0,1/2]Z5-[O1:f:dsp]B3u(a)': 1.0

        #'[0,0,0]GM5+[Ca1:h:dsp]E(a)':1.0,
        #'[0,0,0]GM5+[Ca1:h:dsp]E(b)':1.0,
        #'[0,0,0]GM5-[Cs1:d:dsp]Eu(a)':1.0,
        #'[0,0,0]GM5-[Cs1:d:dsp]Eu(b)':-1.0,
        #'[0,0,0]GM5-[Nb1:a:dsp]Eu(a)':1.0,
        #'[0,0,0]GM5-[Nb1:a:dsp]Eu(b)':1.0,
        #'[0,0,0]GM5-[Nb2:g:dsp]E(a)':1.0,

        #'[1/2,1/2,1/2]A1-[O5:i:dsp]B1(a)':1.0

        #'[1/2,1/2,1/2]A2+[O1:f:dsp]B2u(a)':1.0,
        #'[1/2,1/2,1/2]A2+[O5:i:dsp]B2(a)':1.0,

        #'[1/2,1/2,1/2]A3+[O1:f:dsp]B3u(a)':1.0,
        #'[1/2,1/2,1/2]A3+[O5:i:dsp]B1(a)':1.0,

        #'[1/2,1/2,1/2]A5+[O1:f:dsp]B1u(a)':1.0,
        '[1/2,1/2,1/2]A5+[O1:f:dsp]B1u(b)': 1.0,

        #'[1/2,1/2,0]M1-[O5:i:dsp]B1(a)':1.0,
        #'[1/2,1/2,0]M3+[O1:f:dsp]B3u(a)': 1.0,
        #'[1/2,1/2,0]M3+[O5:i:dsp]B1(a)':1.0,

        #'[1/2,1/2,0]M5+[O1:f:dsp]B1u(a)':1.0,

        #'[0,1/2,0]X3+[O3:g:dsp]E(a)':1.0,
        #'[0,1/2,0]X3+[O3:g:dsp]E(b)':1.0

        #'[0,1/2,1/2]R3+[O1:f:dsp]B1u(a)':1.0,
        #'[0,1/2,1/2]R3+[O1:f:dsp]B1u(b)':1.0,
        #'[0,1/2,1/2]R3+[O3:g:dsp]E(b)': 1.0,
        #'[0,1/2,1/2]R3+[O5:i:dsp]A1(b)':1.0,
    })


def test_parser():
    myparser = Isomode('./isodistort_modes.txt')
    #print(myparser.get_all_mode_fullname())
    print(myparser.symbols)
    myparser.substitute_element({'Cs':'K'})
    #print(myparser.get_mode_displacement("[0,0,1/2]Z5-[O9:g:dsp]E(a)"))
    myparser.read_mode_amplitudes('./mtable.org')
    myparser.prepare_structure_single(
        primitive=True,
        cif_dir='single_mode',
        json_file_name='single_mode.json')


test_parser()
