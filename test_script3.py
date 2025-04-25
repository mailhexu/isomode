"""Test ISODISTORT step by step"""
import os
import re
import tempfile
import requests
from bs4 import BeautifulSoup as BS
from ase import Atoms
import numpy as np
import logging
import time
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

def test_isodistort_step():
    # Create test directory
    TEST_DIR = tempfile.mkdtemp()
    print(f"Using test directory: {TEST_DIR}")

    # Create a simple test structure - cubic perovskite
    # SrTiO3 structure: Sr at corners, Ti at center, O at face centers
    a = 3.9  # lattice parameter in Å
    atoms = Atoms('SrTiO3',
                positions=[[0, 0, 0],         # Sr at origin
                          [a/2, a/2, a/2],    # Ti at center
                          [a/2, a/2, 0],      # O at face centers
                          [a/2, 0, a/2],
                          [0, a/2, a/2]],
                cell=[a, a, a],
                pbc=True)

    # Save structures
    print("\nSaving structures...")
    parent_cif = os.path.join(TEST_DIR, 'parent.cif')
    atoms.write(parent_cif)
    
    distorted = atoms.copy()
    distorted.positions[2,0] += 0.1  # Move O atom along x
    distorted_cif = os.path.join(TEST_DIR, 'distorted.cif')
    distorted.write(distorted_cif)

    # Create session with proper settings
    session = requests.Session()
    session.verify = False
    session.headers.update({
        'User-Agent': 'Mozilla/5.0',
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
        'Connection': 'keep-alive'
    })

    # Step 1: Upload parent structure
    print("\nStep 1: Uploading parent structure...")
    with open(parent_cif, 'rb') as f:
        files = {'toProcess': f}
        ret = session.post(
            "https://stokes.byu.edu/iso/isodistortuploadfile.php",
            data={'input': 'uploadparentcif'},
            files=files)
    
    # Handle automatic form submission
    print("Initial parent response:", ret.text[:200])
    soup = BS(ret.text, 'lxml')
    form = soup.find('form')
    if not form:
        raise RuntimeError("Could not find form in parent upload response")
        
    # Get form data
    data = {}
    for inp in form.find_all('input'):
        name = inp.get('name')
        value = inp.get('value')
        if name and value:
            data[name] = value
            
    print("Parent form data:", data)
    ret = session.post(
        "https://stokes.byu.edu/iso/isodistortform.php",
        data=data)
    
    soup = BS(ret.text, 'lxml')
    title = soup.find('title')
    if title:
        print(f"Parent page title: {title.string}")
    time.sleep(2)

    # Step 2: Upload distorted structure
    print("\nStep 2: Uploading distorted structure...")
    with open(distorted_cif, 'rb') as f:
        files = {'toProcess': f}
        # Add input parameter to initial upload
        data = {'input': 'uploadsubgroupcif'}
        ret = session.post(
            "https://stokes.byu.edu/iso/isodistortuploadfile.php",
            data=data,
            files=files)
    
    # Extract and process the form
    soup = BS(ret.text, 'lxml')
    form = soup.find('form')
    if not form:
        print("Initial response for distorted upload:", ret.text[:200])
        raise RuntimeError("Could not find form in distorted upload response")
        
    # Get the form action URL and data
    action = form.get('action', 'isodistortform.php')
    if not action.startswith('http'):
        action = "https://stokes.byu.edu/iso/" + action
        
    data = {}
    for inp in form.find_all('input'):
        name = inp.get('name')
        value = inp.get('value')
        if name and value:
            data[name] = value
            
    # Add input parameter for subgroup upload
    data['input'] = 'uploadsubgroupcif'
    print("First distorted form data:", data)
    
    ret = session.post(action, data=data)
    
    print("Initial distorted response:")
    print(ret.text[:500])
    
    # Look for basis selection 
    soup = BS(ret.text, 'lxml')
    select = soup.find('select')
    if select:
        print("\nFound basis selection dropdown!")
        options = select.find_all('option')
        print("Available basis options:", [o.get('value') for o in options])
        
        # Get parent form with all data
        basis_form = select.find_parent('form')
        if basis_form:
            data = {}
            for inp in basis_form.find_all('input'):
                name = inp.get('name')
                value = inp.get('value')
                if name and value and name != 'None':
                    data[name] = value
            
            # Add basis selection parameters
            data.update({
                "input": "list",
                "inputbasis": "list",
                "chooseorigin": "false",
                "trynearest": "true"
            })
            
            if len(options) > 1:
                data["basisselect"] = options[1].get('value')
            else:
                data["basisselect"] = options[0].get('value')
            
            print("\nSubmitting basis selection with data:", data)
            ret = session.post(
                "https://stokes.byu.edu/iso/isodistortform.php",
                data=data)
            
            # Check for success
            soup = BS(ret.text, 'lxml')
            title = soup.find('title')
            if title:
                print(f"Response title:", soup.find('title').string if soup.find('title') else "No title")
            print("Response preview:", ret.text[:200])
    else:
        print("Could not find basis selection dropdown")
        print("Full response:", ret.text)
        

    # Step 4: Initialize distortion
    print("\nStep 4: Initializing distortion...")
    
    if 'ISODISTORT: distortion' in ret.text:
        soup = BS(ret.text, 'lxml')
        forms = soup.find_all('form')
        print(f"Found {len(forms)} forms")
        
        for form in forms:
            print("\nForm action:", form.get('action'))
            buttons = form.find_all('input', {'type': 'submit'})
            print("Buttons:", [b.get('value') for b in buttons])
            
            # Look for distortion-related form
            if any(word in str(form).lower() 
                   for word in ['distortion', 'mode', 'irrep']):
                data = {}
                for inp in form.find_all('input'):
                    name = inp.get('name')
                    value = inp.get('value')
                    if name and value and name != 'None':
                        data[name] = value
            
                data.update({
                    "input": "distort",
                    "origintype": "method4",
                    "adjustsub": "false",
                    "primecell": "true",
                    "showstrain": "true",
                    "includerot": "false",
                    "includedisplacive001": "true",
                    "includedisplacive002": "true",
                    "includedisplacive003": "true",
                    "includestrain": "true",
                    "modesdetail": "true"
                })
            
                print("Distortion initialization data:", data)
                #ret = session.post(
                #    "https://stokes.byu.edu/iso/isodistortform.php",
                #    data=data)
            
                #soup = BS(ret.text, 'lxml')
                #title = soup.find('title')
                #print(f"Final page title: {title.string if title else 'No title'}")
                #print("Final response preview:", ret.text[:200])
                break
    else:
        print("Did not find ISODISTORT: distortion page")

if __name__ == "__main__":
    test_isodistort_step()
