#!/usr/bin/env python

from ase.io import read, write
from isomode.pydistort import isocif, isodistort
import argparse
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
    iso.get_mode_amplitude_text()
    mode_details = iso.get_mode_details(save_fname=out_fname)
    return mode_details


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--parent", help="parent cif file name")
    parser.add_argument("-d", "--distort", help="distorted cif file name")
    parser.add_argument(
        "-o",
        "--output",
        help="mode details output filename",
        default="mode_detail.txt")
    args = parser.parse_args()

    view_distort(
        parent_fname=args.parent,
        distorted_fname=args.distort,
        out_fname=args.output)

if __name__=='__main__':
    main()
