#!/usr/bin/env python
from isomode.pydistort import isocif
from ase.io import read
from tempfile import NamedTemporaryFile 
import argparse

def gen_sym_cif(infile, outfile):
    atoms = read(infile)
    # save the input structure to a temporary cif file  
    with NamedTemporaryFile(delete=False, suffix='.cif') as temp_file:
        temp_file_path = temp_file.name
        atoms.write(temp_file_path, format='cif')   

    iso = isocif(temp_file_path)    
    iso.upload_cif()
    iso.findsym()
    iso.save_cif(fname=outfile)

def main():
    parser = argparse.ArgumentParser(
        description='Generate symmetrized CIF file using ISODISTORT'
    )
    parser.add_argument(
        'infile',
        help='Input  file path, does not need to be a CIF file. It can be any ASE-supported file format.'
    )
    parser.add_argument(
        'outfile',
        help='Output CIF file path'
    )
    
    args = parser.parse_args()
    gen_sym_cif(args.infile, args.outfile)

if __name__ == '__main__':
    main()
