#!/usr/bin/env python
import re
import requests
from bs4 import BeautifulSoup as BS
from collections import OrderedDict


from ase.io import read, write
try:
    import spglib
except Exception:
    pass
import tempfile
import os


def tocif(fname, outfname):
    atoms = read(fname)
    atoms.set_pbc([True, True, True])
    write(outfname, atoms)


def view_distort(parent_fname, distorted_fname, out_fname):
    # temp directory
    tmpdir = tempfile.mkdtemp()

    # convert file to cif
    parent_cif = os.path.join(tmpdir, 'parent.cif')
    tocif(parent_fname, outfname=parent_cif)

    # convert highsym_fname
    isosym = isocif(parent_cif)
    isosym.upload_cif()
    isosym.findsym()
    parent_sym_cif = os.path.join(tmpdir, 'parent_sym.cif')
    isosym.save_cif(fname=parent_sym_cif)

    distorted_cif = os.path.join(tmpdir, 'distorted.cif')
    tocif(distorted_fname, outfname=distorted_cif)

    iso = isodistort(parent_cif=parent_sym_cif, distorted_cif=distorted_cif)
    ampt = iso.get_mode_amplitude_text()
    #iso.get_mode_amplitude_text()
    mode_details = iso.get_mode_details(save_fname=out_fname)
    return mode_details

def view_spacegroup(filename='POSCAR',symprec=1e-3):
    atoms=read(filename)
    print("%20s: %s"%('SPACEGROUP', spglib.get_spacegroup(atoms,symprec=symprec)))



class isocif(object):
    def __init__(self, fname):
        self.fname=fname
    def upload_cif(self):
        data = {'input': 'uploadcif'}
        files = {'toProcess': open(self.fname, 'rb')}
        ret = requests.post(
            "http://stokes.byu.edu/iso/isocifuploadfile.php",
            data=data,
            files=files,
            allow_redirects=True)
        text = str(ret.text)

        soup = BS(text,'lxml')
        inputs = soup.find_all('input')
        #print(inputs)
        data = {}
        for i in inputs:
            try:
                name = i.get('name')
                value = i.get('value')
                t = i.get('type')
                if t == 'hidden':
                    data[name] = value
            except:
                pass
        ret = requests.post(
            "http://stokes.byu.edu/iso/isocifform.php", data=data)
        text = ret.text
        self.upload_cif_text = text

    def findsym(self):
        soup = BS(self.upload_cif_text, 'lxml')
        inputs = soup.find_all('input')
        data = {}
        for inp in inputs:
            try:
                name = inp.get('name')
                value = inp.get('value')
                t = inp.get('type')
                data[name] = value
            except:
                pass
        # select basis options. eg.
        data["input"] = "findsym"
        ret = requests.post(
            "http://stokes.byu.edu/iso/isocifform.php", data=data)
        text = ret.text
        self.upload_cif_text = text

    def save_cif(self, fname):
        soup = BS(self.upload_cif_text, 'lxml')
        inputs = soup.find_all('input')
        data = {}
        for inp in inputs:
            try:
                name = inp.get('name')
                value = inp.get('value')
                t = inp.get('type')
                data[name] = value
            except:
                pass
        # select basis options. eg.
        data["input"] = "savecif"
        data["nonstandardsetting"]='false'
        ret = requests.post(
            "http://stokes.byu.edu/iso/isocifform.php", data=data)
        text = ret.text
        self.upload_cif_text = text
        if fname is not None:
            with open(fname, 'w') as myfile:
                myfile.write(text)
        return text



class isodistort(object):
    def __init__(
            self,
            parent_cif="prim_sym.cif",
            distorted_cif="A_0.cif"):
        self.parent_cif = parent_cif
        self.distorted_cif = distorted_cif
        self.upload_parent_cif()
        self.upload_distorted_cif()
        self.select_basis()

    def upload_parent_cif(self):
        data = {'input': 'uploadparentcif'}
        files = {'toProcess': open(self.parent_cif, 'rb')}
        ret = requests.post(
            "http://stokes.byu.edu/iso/isodistortuploadfile.php",
            data=data,
            files=files,
            allow_redirects=True)
        text = str(ret.text)

        fname = re.findall(r'/tmp.*isodistort_.*.iso', text)[0]
        data = {'input': 'uploadparentcif', 'filename': fname}
        ret = requests.post(
            "http://stokes.byu.edu/iso/isodistortform.php",
            data=data,
            allow_redirects=True)
        text = ret.text
        self.upload_parent_cif_text = text

    def upload_distorted_cif(self):
        soup = BS(self.upload_parent_cif_text, 'lxml')
        form_method4 = str(soup.find_all('form')[6])
        #print(form_method4)
        #soup4=BS(form_method4)
        inputs = soup.find_all('input')
        #print(inputs)
        data = {}
        for i in inputs:
            try:
                name = i.get('name')
                value = i.get('value')
                t = i.get('type')
                if t == 'hidden':
                    data[name] = value
            except:
                pass
        #print(data)
        files = {'toProcess': open(self.distorted_cif, 'rb')}

        ret = requests.post(
            "http://stokes.byu.edu/iso/isodistortuploadfile.php",
            data=data,
            files=files,
            allow_redirects=True)
        text = ret.text

        soup = BS(text,'lxml')
        inputs = soup.find_all('input')
        #print(inputs)
        data = {}
        for i in inputs:
            try:
                name = i.get('name')
                value = i.get('value')
                t = i.get('type')
                if t == 'hidden':
                    data[name] = value
            except:
                pass
        data['input'] = 'uploadsubgroupcif'
        ret = requests.post(
            "http://stokes.byu.edu/iso/isodistortform.php", data=data)
        text = ret.text
        self.upload_distorted_cif_text = text

    def select_basis(self):
        # select basis .
        soup = BS(self.upload_distorted_cif_text, 'lxml')
        inputs = soup.find_all('input')
        data = {}
        for inp in inputs:
            try:
                name = inp.get('name')
                value = inp.get('value')
                t = inp.get('type')
                data[name] = value
            except:
                pass
        # select basis options. eg.
        options = soup.find_all('option')
        data["inputbasis"] = "list"
        data["basisselect"] = options[1].get('value')
        data["chooseorigin"] = "false"
        data["trynearest"] = "true"
        #print(data)

        ret = requests.post(
            "http://stokes.byu.edu/iso/isodistortform.php", data=data)
        text = ret.text
        #print(text)
        self.select_basis_text = text

    def get_mode_amplitude_text(self):
        text = self.select_basis_text
        lines = text.split('\n')
        inside = False
        amp_lines = []
        for line in lines:
            if (not inside) and line.find("ampfilename") != -1:
                inside = True
            elif inside:
                if line.find(r"Parameters") != -1:
                    inside = False
                else:
                    amp_lines.append(line)
            else:
                pass
        ret = '\n'.join(amp_lines)
        return ret

    def get_mode_details(self, save_fname=None):
        text = self.select_basis_text
        soup = BS(text, "lxml")
        inputs = soup.find_all('input')
        data = {}
        for inp in inputs:
            try:
                name = inp.get('name')
                value = inp.get('value')
                t = inp.get('type')
                if t not in ['radio', 'checkbox']:
                    data[name] = value
            except:
                pass
        # select basis options. eg.
        options = soup.find_all('option')
        data["topasstrain"] = "false"
        data["treetopas"] = "false"
        data["cifmovie"] = "false"
        data["nonstandardsetting"] = 'false'
        data["origintype"] = "modesdetails"
        data["varcifmovie"] = 'linear'
        data["cifdec"] = " 5"
        ret = requests.post(
            "http://stokes.byu.edu/iso/isodistortform.php", data=data)
        text = ret.text
        lines = text.split('\n')
        p = re.compile(r'<pre>([\s|\S]*)<\/pre>', re.MULTILINE)
        texts = p.findall(text)
        if len(texts) != 0:
            text = texts[0]
        else:
            text = ''
        soup=BS(text, "html.parser")
        text=soup.get_text()    
        if save_fname is not None:
            with open(save_fname,'w') as myfile:
                myfile.write(text)
        return text

def get_summary(text, remove0=True):
    modes=OrderedDict()
    for line in text.splitlines():
        if 'all' in line:
            if 'Overall' in line:
                total=(float(line.strip().split()[1]),float(line.strip().split()[2]))
            else:
                s=line.strip().split()
                if remove0==False or abs(float(s[2]))>0.00001:
                     modes[s[0]]=(float(s[2]), float(s[3]))
    return total, modes

def print_summary(text):
    total, modes = get_summary(text)
    for m in sorted(modes.keys(), key=lambda x: modes[x][0], reverse=True):
        v=modes[m]
        print("%20s:  %.4f  %.4f"%(m, v[0],  v[1]))
    print("%20s:  %.4f  %.4f"%('Total',total[0], total[1]))

def test_isocif(fname='nmodes/primitive.cif'):
    iso = isocif(fname)
    iso.upload_cif()
    iso.findsym()
    iso.save_cif(fname='save.cif')

def test(parent_cif='save.cif', distorted_cif='nmodes/A_0.cif', mode_detail_file='mode_detail.txt'):
    iso = isodistort(parent_cif=parent_cif, distorted_cif=distorted_cif)
    ampt = iso.get_mode_amplitude_text()
    mode_details=iso.get_mode_details(save_fname='mode_detail.txt')

#test_isocif()
#test()
