#!/usr/bin/env python3
"""
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

Created on Mon Oct 25 11:50:39 2021

@author: mollygravett
modified jonathan pickering 23Aug22
"""
# set up linting conditions
# pylint: disable = import-error
import argparse
import pathlib
from enum import Enum
import mrcfile
from scipy import ndimage
import numpy as np

class Algorithm(Enum):
    """
    enumeration of possible kernels
    """
    ELLIPSOID = 'ellipsoid'
    GAUSSIAN = 'gaussian'
    UNIFORM = 'uniform'

def get_args():
    """
    get the command line arguments
        Returns
            (argparse.namespace)
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("-i",
                        "--input",
                        type=pathlib.Path,
                        required=True,
                        help="input file")

    parser.add_argument("-o",
                        "--output",
                        type=pathlib.Path,
                        required=True,
                        help="output file")

    parser.add_argument("-x",
                        "--xstart",
                        type=int,
                        default=0,
                        help="start of x section")

    parser.add_argument("-X",
                        "--xend",
                        type=int,
                        default=100,
                        help="end of x section")

    return parser.parse_args()


def output_mrcfile(data, out_file, original_mrc):
    """
    write out a scaled mrc file
        Args:
            data (numpy.ndarray): the volume data
            out_file (pathlib.Path) the output file
            original_mrc (mrcfile.File) the original input file
    """
    with mrcfile.new(out_file, overwrite=True) as mrc:
        mrc.set_data(data)
        mrc.header.origin = original_mrc.header.origin
        mrc.voxel_size = original_mrc.voxel_size

def chop_image(mrc, args):
    print(f"chop from {args.xstart} to {args.xend}")
    pass

def main():
    """
    run the script
    """
    args = get_args()
    with mrcfile.open(args.input, mode='r') as mrc:
        chop_image(mrc, args)
    print(f"{args.input} filtered and written to {args.output}")

if __name__ == "__main__":
    main()
