#!/usr/bin/env python
from isomode.pydistort import isocif
from ase.io import read, write
import sys

def gen_sym_cif(infile, outfile):
    atoms=read(infile)
    iso = isocif(fname)
    iso.upload_cif()
    iso.findsym()
    iso.save_cif(fname=outfile)


if __name__=='__main__':
    infile=sys.argv[1]
    outfile=sys.argv[2]

