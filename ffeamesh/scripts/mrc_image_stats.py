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

 @author J Pickering & J Leng Universary of Leeds 2022
"""
# set up linting conditions
# pylint: disable = import-error

import argparse
import pathlib
import sys
from scipy import stats
import numpy as np
import mrcfile

def get_args():
    """
    get the command line arguments
        Returns
            (argparse.namespace)
    """
    parser = argparse.ArgumentParser("""read a mrc file and print out image intensity stats""")

    parser.add_argument("-i",
                        "--input",
                        type=pathlib.Path,
                        required=True,
                        help="input file")

    parser.add_argument("-o",
                        "--output",
                        type=pathlib.Path,
                        required=False,
                        help="output file")

    parser.add_argument("-w",
                        "--overwrite",
                        action="store_true",
                        help="overwrite")

    return parser.parse_args()

def validate_args(args):
    """
    check the command line arguments are sane
    Args
        args (argparse.Namespace)
    Returns:
        (str/None) error message if problem else None
    """
    if not args.input.exists():
        return f"file {args.input} does not exist!"

    if args.output is not None and args.output.exists():
        if not args.overwrite:
            return f"file {args.output} already exists run with -w to overwrite"

    return None

def print_image_stats(mrc, infile, outfile=None):
    """
    generate and print the stats of the image array
    Args:
        mrc (mrcfile): the image source
        infile (str): file name of input
        outfile (pathlib.Path) the output file
    """

    descriptive_stats = stats.describe(mrc.data.flatten())
    histogram = np.histogram(mrc.data.flatten(), bins=10)
    size = mrc.header.mx * mrc.header.my * mrc.header.mz

    print(f"statistics for file {infile}")
    print(f"Number of voxels in image {size}")
    print(descriptive_stats)
    print(f"bins   {histogram[1]}")
    print(f"counts {histogram[0]}")

def main():
    """
    run the script
    """
    args = get_args()
    message = validate_args(args)
    if message is not None:
        print(message, file=sys.stderr, flush=True)
        sys.exit()

    with mrcfile.mmap(args.input, 'r') as mrc:
        print_image_stats(mrc, str(args.input), args.output)

if __name__ == "__main__":
    main()
