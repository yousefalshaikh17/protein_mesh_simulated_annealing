#!/usr/bin/env python
"""
 mrc_resample.py

 A script to coarsen an MRC image file from 5Å to 15Å for example.

 Replaces the zoom function developed by Molly Gravett.

 -------------------------------------

 This file is part of the FFEA simulation package

 Copyright (c) by the Theory and Development FFEA teams,
 as they appear in the README.md file.

 FFEA is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 FFEA is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with FFEA.  If not, see <http://www.gnu.org/licenses/>.

 To help us fund FFEA development, we humbly ask that you cite
 the research papers on the package.

 Created on Wed 01-Nov-23

 @author: Jonathan Pickering
"""
# set up linting conditions
# pylint: disable = import-error

import argparse
import pathlib

from ffeamesh.mrc_zoom import resample_volume_data

def get_args():
    """
    get the command line arguments
    Returns
        (argparse.namespace)
    """
    parser = argparse.ArgumentParser("""coarsen an MRC image file.""" )

    parser.add_argument("-i",
                        "--input",
                        type=existing_mrc_file,
                        required=True,
                        help="input file")

    parser.add_argument("-o",
                        "--output",
                        type=mrc_file,
                        required=True,
                        help="output file")

    parser.add_argument("-s",
                        "--scalefactor",
                        type=valid_scalefactor,
                        default=0.5,
                        help="scalefactor by which to coarsen input MRC.")

    return parser.parse_args()

def mrc_file(filename):
    """
    make a pathlib.Path and check it has postfix .mrc
    Args:
        filename (str): name/path of file
    Returns:
        pathlib.Path
    Raises:
        ValueError if the filename does not end .mrc
    """
    path = pathlib.Path(filename)
    return path.with_suffix('.mrc')

def existing_mrc_file(filename):
    """
    make a pathlib.Path and check it exists
    Args:
        filename (str): name/path of file
    Returns:
        pathlib.Path
    Raises:
        ValueError if filename cannot be made into path
    """
    path = mrc_file(filename)

    if not path.exists():
        raise ValueError(f"File {filename} does not exist!")

    return path

def valid_scalefactor(scalefactor_s):
    """
    ensure the scalefactor is sane
    Args:
        scalefactor_s (str): scalefactor in string form
    Returns
        (0.0<float<1.0): the scalefactor
    Raises:
        ValueError if can not convert to float or value <= 0.0
    """
    scalefactor = float(scalefactor_s)

    if scalefactor <= 0.0 or scalefactor >= 1.0:
        raise ValueError(f"Error {scalefactor} outside range (0.0, 1.0)")

    return scalefactor

def main():
    """
    run the script
    """
    args = get_args()
    # TODO add warning if very small scale

    resample_volume_data(args.input, args.output, args.scalefactor)

    # TODO output the numbers & size of voxels
    print(f"{args.input} scaled to {args.scalefactor} and written to {args.output}")

if __name__ == "__main__":
    main()
