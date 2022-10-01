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
        return f"input file {args.input} does not exist!"

    if args.output is not None and args.output.exists():
        if not args.overwrite:
            return f"output file {args.output} already exists run with -w to overwrite"

    return None

def get_image_stats(mrc):
    """
    generate the stats of the image array
    Args:
        mrc (mrcfile): the image source
    Returns:
        (scipy.stats.DescribeResult)
        (numpy.histogram)
    """
    descriptive_stats = stats.describe(mrc.data.flatten())
    histogram = np.histogram(mrc.data.flatten(), bins=10)

    return descriptive_stats, histogram

def print_image_stats(descriptive_stats,  histogram, infile, outfile=None):
    """
    generate and print the stats of the image array
    Args:
        descriptive_stats (scipy.stats.DescribeResult): mean etc
        histogram (numpy.histogram): ten bin histogram
        mrc (mrcfile): the image source
        infile (str): file name of input
        outfile (pathlib.Path) the output file
    """

    print(f"File, {infile}", file=outfile)
    print(f"Number of voxels in image, {descriptive_stats.nobs}", file=outfile)
    print(f"Minimum, {descriptive_stats.minmax[0]:.6}", file=outfile)
    print(f"Maximum, {descriptive_stats.minmax[1]:.6}", file=outfile)
    print(f"Mean, {descriptive_stats.mean:.6}", file=outfile)
    print(f"Variance {descriptive_stats.variance:.6}", file=outfile)

    print("\nbin, range, count", file=outfile)
    for i, count in enumerate(histogram[0]):
        print(f"{i}, ({histogram[1][i]:.6}, {histogram[1][i+1]:.6}), {count}", file=outfile)

def main():
    """
    run the script
    """
    args = get_args()
    message = validate_args(args)
    if message is not None:
        print(f"Error: {message}", file=sys.stderr, flush=True)
        sys.exit()

    with mrcfile.mmap(args.input, 'r') as mrc:
        descriptive_stats,  histogram = get_image_stats(mrc)
        if args.output is not None:
            with args.output.open("w") as outfile:
                print_image_stats(descriptive_stats,  histogram, str(args.input), outfile)
        else:
            print_image_stats(descriptive_stats,  histogram, str(args.input))

if __name__ == "__main__":
    main()
