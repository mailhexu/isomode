#!/usr/bin/env python
import argparse
from isomode.pydistort import isocif, isodistort, print_summary, view_distort, view_spacegroup

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


    #raise NotImplementedError(
    #    "This script is not working. "
    #)   

    mode_details=view_distort(
        parent_fname=args.parent,
        distorted_fname=args.distort,
        out_fname=args.output)
    
    print_summary(mode_details)

if __name__=='__main__':
    main()
