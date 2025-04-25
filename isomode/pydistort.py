#!/usr/bin/env python
"""
view the amplitude of distortions with isotropy
"""
import re
import requests
from bs4 import BeautifulSoup as BS
from collections import OrderedDict
import logging
from ase.io import read, write
import urllib.parse
try:
    import spglib
except Exception:
    pass
import tempfile
import os
import time
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

def tocif(fname, outfname):
    atoms = read(fname)
    atoms.set_pbc([True, True, True])
    write(outfname, atoms)

class isocif(object):
    def __init__(self, fname):
        self.fname=fname
        self.session = requests.Session()
        self.session.verify = False
        self.session.headers.update({
            'User-Agent': 'Mozilla/5.0',
            'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8'
        })

    def upload_cif(self):
        logger.info("Uploading CIF...")
        data = {'input': 'uploadcif'}
        with open(self.fname, 'rb') as f:
            files = {'toProcess': f}
            ret = self.session.post(
                "https://stokes.byu.edu/iso/isocifuploadfile.php",
                data=data,
                files=files)
        text = str(ret.text)
        soup = BS(text, 'lxml')
        inputs = soup.find_all('input')
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
        ret = self.session.post(
            "https://stokes.byu.edu/iso/isocifform.php", data=data)
        text = ret.text
        self.upload_cif_text = text

    def findsym(self):
        logger.info("Finding symmetry...")
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
        data["input"] = "findsym"
        ret = self.session.post(
            "https://stokes.byu.edu/iso/isocifform.php", data=data)
        text = ret.text
        self.upload_cif_text = text
        
        soup = BS(text, 'lxml')
        spacegroup_text = soup.find(text=re.compile(r'Space Group:.*'))
        if spacegroup_text:
            logger.info(f"Detected {spacegroup_text.strip()}")

    def save_cif(self, fname):
        logger.info(f"Saving CIF to {fname}")
        soup = BS(self.upload_cif_text, 'lxml')
        inputs = soup.find_all('input')
        data = {}
        for inp in inputs:
            try:
                name = inp.get('name')
                value = inp.get('value')
                t = inp.get('type')
                if t == 'hidden':
                    data[name] = value
            except:
                pass
        data["input"] = "savecif"
        data["nonstandardsetting"]='false'
        ret = self.session.post(
            "https://stokes.byu.edu/iso/isocifform.php", data=data)
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
            distorted_cif=None):
        self.parent_cif = parent_cif
        self.distorted_cif = distorted_cif
        self.session = requests.Session()
        self.session.verify = False
        self.session.headers.update({
            'User-Agent': 'Mozilla/5.0',
            'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
            'Connection': 'keep-alive'
        })
        logger.info("Initializing ISODISTORT analysis")
        # Upload parent structure and process it
        self.upload_parent_cif()
        if distorted_cif:
            # Upload distorted structure and process it
            self.upload_distorted_cif()
            # Select basis
            self.select_basis()
            # Initialize distortion analysis to get mode details
            self.initialize_distortion()

    def upload_parent_cif(self):
        logger.info("Uploading parent structure")
        data = {'input': 'uploadparentcif'}
        with open(self.parent_cif, 'rb') as f:
            files = {'toProcess': f}
            ret = self.session.post(
                "https://stokes.byu.edu/iso/isodistortuploadfile.php",
                data=data,
                files=files)
            
        text = str(ret.text)
        matches = re.findall(r'/tmp.*isodistort_.*.iso', text)
        if not matches:
            raise RuntimeError("Failed to get temporary file name")
            
        fname = matches[0]
        data = {'input': 'uploadparentcif', 'filename': fname}
        
        ret = self.session.post(
            "https://stokes.byu.edu/iso/isodistortform.php",
            data=data)
        text = ret.text
        self.upload_parent_cif_text = text

    def upload_distorted_cif(self):
        logger.info("Uploading distorted structure")
        soup = BS(self.upload_parent_cif_text, 'lxml')
        inputs = soup.find_all('input')
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
        files = {'toProcess': open(self.distorted_cif, 'rb')}
        ret = self.session.post(
            "https://stokes.byu.edu/iso/isodistortuploadfile.php",
            data=data,
            files=files)
        text = ret.text

        soup = BS(text,'lxml')
        inputs = soup.find_all('input')
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
        ret = self.session.post(
            "https://stokes.byu.edu/iso/isodistortform.php", data=data)
        text = ret.text
        self.upload_distorted_cif_text = text

    def select_basis(self):
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

        options = soup.find_all('option')
        data["inputbasis"] = "list"
        data["basisselect"] = options[1].get('value')
        data["chooseorigin"] = "false"
        data["trynearest"] = "true"
        print("Basis selection data:", data)

        ret = self.session.post(
            "https://stokes.byu.edu/iso/isodistortform.php", data=data)
        text = ret.text
        self.select_basis_text = text
        print("--- Raw response from basis selection after select basis ---")
        print(text)
        print("----------------------------------------")

    def get_distortion_details(self):
        """Get distortion details after basis selection, this is by checking the "mode details: option and click the OK button"""
        logger.info("Getting distortion details")
        soup = BS(self.select_basis_text, 'lxml')
        form = soup.find('form', action='isodistortform.php')
        if not form:
            print("Initial response for distortion details:", self.select_basis_text[:200])
            raise RuntimeError("Could not find form in distortion details response")
        data = {}
        if form:
            inputs = form.find_all('input', type='hidden')
            for inp in inputs:
                name = inp.get('name')
                value = inp.get('value')
                #if name == 'mapatomsdata':
                #    print("--- Raw mapatomsdata string ---")
                #    print(value)
                #    print("-----------------------------")
                #    # Encode mapatomsdata to Latin-1
                #    data[name] = value.encode('latin-1')
                #else:
                #    data[name] = value

        data["origintype"] = "method4"


        headers = self.session.headers.copy()
        headers['Content-Type'] = 'application/x-www-form-urlencoded'

        ret = self.session.post(
            "https://stokes.byu.edu/iso/isodistortform.php", data=data, headers=headers)
        text = ret.text

    def initialize_distortion(self):
        """Initialize distortion analysis"""
        logger.info("Initializing distortion analysis")
        soup = BS(self.select_basis_text, 'lxml')
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

        data["topasstrain"] = "false"
        data["treetopas"] = "false"
        data["cifmovie"] = "false"
        data["nonstandardsetting"] = 'false'
        data["origintype"] = "modesdetails"
        data["varcifmovie"] = 'linear'
        data["cifdec"] = " 5"
        ret = self.session.post(
            "https://stokes.byu.edu/iso/isodistortform.php", data=data)
        text = ret.text
        print("--- Raw response from first POST in initialize_distortion ---")
        print(text)
        print("-------------------------------------------------------------")

        # Parse the response to find the form for mode details
        soup = BS(text, "lxml")
        form = soup.find('form', action='isodistortform.php')
        
        # use default value for the form except for the 'modesdetails' input
        if form:
            data_for_modes = {}
            inputs = form.find_all('input', type='hidden')
            for inp in inputs:
                name = inp.get('name')
                value = inp.get('value')
                data_for_modes[name] = value

            # Add the 'modesdetails' input as per step 5
            data_for_modes["origintype"] = "modesdetails"

            print("--- Data being sent for mode details ---")
            print(data_for_modes)
            print("----------------------------------------")

            # Make the second POST request to get mode details
            ret_modes = self.session.post(
                "https://stokes.byu.edu/iso/isodistortform.php", data=data_for_modes)
            text_modes = ret_modes.text

            print("--- Raw response from second POST (mode details) ---")
            print(text_modes)
            print("----------------------------------------------------")

            # Extract text from the <pre> tags in the second response
            soup_modes = BS(text_modes, "html.parser")
            pre_tag = soup_modes.find('pre')
            if pre_tag:
                self.mode_details_text = pre_tag.get_text()
            else:
                self.mode_details_text = ""
        else:
            logger.error("Could not find the form to get mode details.")
            self.mode_details_text = ""

    def get_mode_details(self, save_fname=None):
        if not hasattr(self, 'mode_details_text'):
            return ''
        text = self.mode_details_text
        if save_fname is not None:
            with open(save_fname,'w') as myfile:
                myfile.write(text)
        return text

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

def test_isocif(fname='nmodes/primitive.cif'):
    iso = isocif(fname)
    iso.upload_cif()
    iso.findsym()
    iso.save_cif(fname='save.cif')

def test(parent_cif='save.cif', distorted_cif='nmodes/A_0.cif', mode_detail_file='mode_detail.txt'):
    iso = isodistort(parent_cif=parent_cif, distorted_cif=distorted_cif)
    ampt = iso.get_mode_amplitude_text()
    mode_details = iso.get_mode_details(save_fname=mode_detail_file)
    return mode_details
