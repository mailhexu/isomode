#!/usr/bin/env python
"""
view the amplitude of distortions with isotropy
"""
import re
import requests
import logging
import sys
import os
import tempfile
from bs4 import BeautifulSoup as BS
from collections import OrderedDict
from pathlib import Path

from ase.io import read, write
try:
    import spglib
except Exception:
    pass

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler('pydistort.log')
    ]
)
logger = logging.getLogger(__name__)

class PyDistortError(Exception):
    """Base exception class for pydistort errors"""
    pass

class NetworkError(PyDistortError):
    """Raised when network operations fail"""
    pass

class FileOperationError(PyDistortError):
    """Raised when file operations fail"""
    pass

def tocif(fname, outfname):
    """Convert input file to CIF format"""
    try:
        logger.info(f"Converting {fname} to CIF format: {outfname}")
        atoms = read(fname)
        atoms.set_pbc([True, True, True])
        write(outfname, atoms)
        logger.info(f"Successfully converted {fname} to CIF")
    except Exception as e:
        logger.error(f"Failed to convert {fname} to CIF: {str(e)}")
        raise FileOperationError(f"Failed to convert {fname} to CIF: {str(e)}")

def view_distort(parent_fname, distorted_fname, out_fname):
    """View the distortion between parent and distorted structures"""
    try:
        # Create temp directory
        tmpdir = tempfile.mkdtemp()
        logger.info(f"Created temporary directory: {tmpdir}")

        # Convert files to CIF
        logger.info("Converting parent file to CIF")
        parent_cif = os.path.join(tmpdir, 'parent.cif')
        tocif(parent_fname, outfname=parent_cif)

        logger.info("Converting parent CIF file to high symmetry using findsym")
        isosym = isocif(parent_cif)
        isosym.upload_cif()
        isosym.findsym()
        parent_sym_cif = os.path.join(tmpdir, 'parent_sym.cif')
        isosym.save_cif(fname=parent_sym_cif)

        logger.info("Converting distorted file to CIF")
        distorted_cif = os.path.join(tmpdir, 'distorted.cif')
        tocif(distorted_fname, outfname=distorted_cif)


        atoms = read(distorted_cif) 

        print(atoms.get_scaled_positions()) 

        logger.info("Getting distortion mode details")
        iso = isodistort(parent_cif=parent_sym_cif, distorted_cif=distorted_cif)
        ampt = iso.get_mode_amplitude_text()
        mode_details = iso.get_mode_details(save_fname=out_fname)
        return mode_details
    except Exception as e:
        logger.error(f"Error in view_distort: {str(e)}")
        raise PyDistortError(f"Failed to analyze distortion: {str(e)}")

def view_spacegroup(filename='POSCAR', symprec=1e-3):
    """View the space group of a structure"""
    try:
        logger.info(f"Analyzing space group for {filename}")
        atoms = read(filename)
        spacegroup = spglib.get_spacegroup(atoms, symprec=symprec)
        logger.info(f"Space group determined: {spacegroup}")
        print("%20s: %s" % ('SPACEGROUP', spacegroup))
    except Exception as e:
        logger.error(f"Failed to determine space group: {str(e)}")
        raise PyDistortError(f"Failed to determine space group: {str(e)}")

class isocif(object):
    def __init__(self, fname):
        self.fname = fname
        logger.info(f"Initialized isocif with file: {fname}")

    def upload_cif(self):
        """Upload CIF file to server"""
        try:
            logger.info(f"Uploading CIF file: {self.fname}")
            data = {'input': 'uploadcif'}
            with open(self.fname, 'rb') as f:
                files = {'toProcess': (self.fname, f)}  
                ret = requests.post(
                    "https://stokes.byu.edu/iso/isocifuploadfile.php",
                    data=data,
                    files=files,
                    #allow_redirects=True
                    )
            
            text = str(ret.text)
            logger.debug(f"Upload response: {text}")

            soup = BS(text, 'lxml')
            inputs = soup.find_all('input')
            data = {i.get('name'): i.get('value') 
                   for i in inputs 
                   if i.get('type') == 'hidden' and i.get('name')}

            ret = requests.post(
                "https://stokes.byu.edu/iso/isocifform.php", data=data)
            self.upload_cif_text = ret.text
            logger.info("CIF file uploaded successfully")
        except requests.RequestException as e:
            logger.error(f"Network error during CIF upload: {str(e)}")
            raise NetworkError(f"Failed to upload CIF: {str(e)}")
        except Exception as e:
            logger.error(f"Error during CIF upload: {str(e)}")
            raise PyDistortError(f"Failed to upload CIF: {str(e)}")

    def findsym(self):
        """Find symmetry using uploaded CIF"""
        try:
            logger.info("Finding symmetry")
            soup = BS(self.upload_cif_text, 'lxml')
            inputs = soup.find_all('input')
            data = {inp.get('name'): inp.get('value') 
                   for inp in inputs 
                   if inp.get('name')}
            data["input"] = "findsym"
            
            ret = requests.post(
                "https://stokes.byu.edu/iso/isocifform.php", data=data)
            self.upload_cif_text = ret.text
            logger.info("Symmetry analysis completed")
        except requests.RequestException as e:
            logger.error(f"Network error during symmetry analysis: {str(e)}")
            raise NetworkError(f"Failed to find symmetry: {str(e)}")
        except Exception as e:
            logger.error(f"Error during symmetry analysis: {str(e)}")
            raise PyDistortError(f"Failed to find symmetry: {str(e)}")

    def save_cif(self, fname):
        """Save processed CIF to file"""
        try:
            logger.info(f"Saving CIF to {fname}")
            soup = BS(self.upload_cif_text, 'lxml')
            inputs = soup.find_all('input')
            data = {inp.get('name'): inp.get('value') 
                   for inp in inputs 
                   if inp.get('name')}
            data["input"] = "savecif"
            data["nonstandardsetting"] = 'false'
            
            ret = requests.post(
                "https://stokes.byu.edu/iso/isocifform.php", data=data)
            self.upload_cif_text = ret.text
            
            if fname is not None:
                with open(fname, 'w') as myfile:
                    myfile.write(ret.text)
                logger.info(f"CIF file saved to {fname}")
            return ret.text
        except requests.RequestException as e:
            logger.error(f"Network error while saving CIF: {str(e)}")
            raise NetworkError(f"Failed to save CIF: {str(e)}")
        except Exception as e:
            logger.error(f"Error while saving CIF: {str(e)}")
            raise PyDistortError(f"Failed to save CIF: {str(e)}")

class isodistort(object):
    def __init__(
            self,
            parent_cif="prim_sym.cif",
            distorted_cif="A_0.cif"):
        self.parent_cif = parent_cif
        self.distorted_cif = distorted_cif
        logger.info(f"Initialized isodistort with parent: {parent_cif}, distorted: {distorted_cif}")
        self.upload_parent_cif()
        self.upload_distorted_cif()
        self.select_basis()

    def upload_parent_cif(self):
        """Upload parent CIF file"""
        try:
            logger.info(f"Uploading parent CIF: {self.parent_cif}")
            data = {'input': 'uploadparentcif'}
            with open(self.parent_cif, 'rb') as f:
                files = {'toProcess': (self.parent_cif, f)} 
                ret = requests.post(
                    "https://stokes.byu.edu/iso/isodistortuploadfile.php",
                    data=data,
                    files=files,
                    #allow_redirects=True
                    )
            text = str(ret.text)
            logger.debug(f"Upload response: {text}")

            fname = re.findall(r'/tmp.*isodistort_.*.iso', text)[0]
            data = {'input': 'uploadparentcif', 'filename': fname}
            ret = requests.post(
                "https://stokes.byu.edu/iso/isodistortform.php",
                data=data,
                #allow_redirects=True
                )
            self.upload_parent_cif_text = ret.text
            logger.info("Parent CIF uploaded successfully")
        except requests.RequestException as e:
            logger.error(f"Network error uploading parent CIF: {str(e)}")
            raise NetworkError(f"Failed to upload parent CIF: {str(e)}")
        except Exception as e:
            logger.error(f"Error uploading parent CIF: {str(e)}")
            raise PyDistortError(f"Failed to upload parent CIF: {str(e)}")

    def upload_distorted_cif(self):
        """Upload distorted CIF file"""
        try:
            logger.info(f"Uploading distorted CIF: {self.distorted_cif}")
            soup = BS(self.upload_parent_cif_text, 'lxml')
            inputs = soup.find_all('input')
            data = {i.get('name'): i.get('value') 
                   for i in inputs 
                    #if i.get('name')}
                    if i.get('type') == 'hidden' and i.get('name')}

            with open(self.distorted_cif, 'rb') as f:
                files = {'toProcess': (self.distorted_cif, f)}
                ret = requests.post(
                    "https://stokes.byu.edu/iso/isodistortuploadfile.php",
                    data=data,
                    files=files,
                    #allow_redirects=True
                    )
            text = ret.text

            soup = BS(text, 'lxml')
            inputs = soup.find_all('input')
            data = {i.get('name'): i.get('value') 
                   for i in inputs 
                    #if i.get('type') == 'hidden' and i.get('name')}
                    if i.get('name')}
            data['input'] = 'uploadsubgroupcif'

            ret = requests.post(
                "https://stokes.byu.edu/iso/isodistortform.php", data=data, allow_redirects=True)
            self.upload_distorted_cif_text = ret.text
            logger.info("Distorted CIF uploaded successfully")
        except requests.RequestException as e:
            logger.error(f"Network error uploading distorted CIF: {str(e)}")
            raise NetworkError(f"Failed to upload distorted CIF: {str(e)}")
        except Exception as e:
            logger.error(f"Error uploading distorted CIF: {str(e)}")
            raise PyDistortError(f"Failed to upload distorted CIF: {str(e)}")

    def select_basis(self):
        """Select basis for distortion analysis"""
        try:
            logger.info("Selecting basis for analysis")
            soup = BS(self.upload_distorted_cif_text, 'lxml')
            inputs = soup.find_all('input')
            data = {inp.get('name'): inp.get('value') 
                   for inp in inputs 
                   if inp.get('name')}
            
            options = soup.find_all('option')
            data["inputbasis"] = "list"
            data["basisselect"] = options[1].get('value')
            data["chooseorigin"] = False
            data["trynearest"] = True
            data["domapatoms"] = 0
            data["zeromodes"] =  False
            data['input'] = 'distort'
            data['origintype'] = 'method4'

            ret = requests.post(
                "https://stokes.byu.edu/iso/isodistortform.php", data=data, allow_redirects=True)
            self.select_basis_text = ret.text
            print(self.select_basis_text)   
            logger.info("Basis selection completed")
        except requests.RequestException as e:
            logger.error(f"Network error during basis selection: {str(e)}")
            raise NetworkError(f"Failed to select basis: {str(e)}")
        except Exception as e:
            logger.error(f"Error during basis selection: {str(e)}")
            raise PyDistortError(f"Failed to select basis: {str(e)}")

    def get_mode_amplitude_text(self):
        """Extract mode amplitude text from analysis"""
        try:
            logger.info("Extracting mode amplitudes")
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
            ret = '\n'.join(amp_lines)
            logger.info("Mode amplitudes extracted successfully")
            return ret
        except Exception as e:
            logger.error(f"Error extracting mode amplitudes: {str(e)}")
            raise PyDistortError(f"Failed to extract mode amplitudes: {str(e)}")

    def get_mode_details(self, save_fname=None):
        """Get detailed mode analysis"""
        try:
            logger.info("Getting mode details")
            text = self.select_basis_text
            soup = BS(text, "lxml")
            inputs = soup.find_all('input')
            data = {inp.get('name'): inp.get('value')
                   for inp in inputs
                   #if inp.get('type') not in ['radio', 'checkbox'] and inp.get('name')}
                   if inp.get('name')}

            #options = soup.find_all('option')
            data.update({
                "topasstrain": "false",
                "treetopas": "false",
                "cifmovie": "false",
                "nonstandardsetting": 'false',
                "origintype": "modesdetails",
                "varcifmovie": 'linear',
                "zeromodes": 'false',   
                "cifdec": " 5"
            })
            ret = requests.post(
                "https://stokes.byu.edu/iso/isodistortform.php", data=data, 
                #allow_redirects=True
                )
            text = ret.text

            
            p = re.compile(r'<pre>([\s|\S]*)<\/pre>', re.MULTILINE)
            texts = p.findall(text)
            if texts:
                text = texts[0]
            else:
                text = ''
            
            soup = BS(text, "html.parser")
            text = soup.get_text()

            if save_fname is not None:
                with open(save_fname, 'w') as myfile:
                    myfile.write(text)
                logger.info(f"Mode details saved to {save_fname}")
            return text
        except requests.RequestException as e:
            logger.error(f"Network error getting mode details: {str(e)}")
            raise NetworkError(f"Failed to get mode details: {str(e)}")
        except Exception as e:
            logger.error(f"Error getting mode details: {str(e)}")
            raise PyDistortError(f"Failed to get mode details: {str(e)}")

def get_summary(text, remove0=True):
    """Get summary of mode analysis"""
    try:
        logger.info("Generating mode summary")
        modes = OrderedDict()
        for line in text.splitlines():
            if 'all' in line:
                if 'Overall' in line:
                    total = (float(line.strip().split()[1]), float(line.strip().split()[2]))
                else:
                    s = line.strip().split()
                    if remove0 is False or abs(float(s[2])) > 0.00001:
                        modes[s[0]] = (float(s[2]), float(s[3]))
        logger.info("Mode summary generated successfully")
        return total, modes
    except Exception as e:
        logger.error(f"Error generating mode summary: {str(e)}")
        raise PyDistortError(f"Failed to generate mode summary: {str(e)}")

def print_summary(text):
    """Print summary of mode analysis"""
    try:
        logger.info("Printing mode summary")
        total, modes = get_summary(text)
        for m in sorted(modes.keys(), key=lambda x: modes[x][0], reverse=True):
            v = modes[m]
            print("%20s:  %.4f  %.4f" % (m, v[0], v[1]))
        print("%20s:  %.4f  %.4f" % ('Total', total[0], total[1]))
        logger.info("Mode summary printed successfully")
    except Exception as e:
        logger.error(f"Error printing mode summary: {str(e)}")
        raise PyDistortError(f"Failed to print mode summary: {str(e)}")

def test_isocif(fname='nmodes/primitive.cif'):
    """Test isocif functionality"""
    try:
        logger.info(f"Testing isocif with file: {fname}")
        iso = isocif(fname)
        iso.upload_cif()
        iso.findsym()
        iso.save_cif(fname='save.cif')
        logger.info("isocif test completed successfully")
    except Exception as e:
        logger.error(f"Error in isocif test: {str(e)}")
        raise PyDistortError(f"isocif test failed: {str(e)}")

def test(parent_cif='save.cif', distorted_cif='nmodes/A_0.cif', mode_detail_file='mode_detail.txt'):
    """Test full distortion analysis"""
    try:
        logger.info("Running full distortion analysis test")
        iso = isodistort(parent_cif=parent_cif, distorted_cif=distorted_cif)
        ampt = iso.get_mode_amplitude_text()
        mode_details = iso.get_mode_details(save_fname=mode_detail_file)
        logger.info("Distortion analysis test completed successfully")
    except Exception as e:
        logger.error(f"Error in distortion analysis test: {str(e)}")
        raise PyDistortError(f"Distortion analysis test failed: {str(e)}")

if __name__ == "__main__":
    try:
        test_isocif()
        test()
    except Exception as e:
        logger.error(f"Error in main execution: {str(e)}")
        sys.exit(1)
