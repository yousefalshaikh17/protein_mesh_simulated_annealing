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

 DESCRIPTION
 to coarsen a mesh (to 15 Å for example) in chimerax use:
 vol resample #1 spacing 15 - this coarsens to 15 Å voxels
 save newmap.mrc model #2 - this saves your new coarsened model as "newmap.mrc"

 @author: mollygravett
 modified jonathan pickering 23Aug22
"""



import argparse
import pathlib
import sys

from ffeamesh.mrc_zoom import refine_volume_data

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

    parser.add_argument("-r",
                        "--resolution",
                        type=float,
                        default=20.0,
                        help="Resolution by which to coarsen input MRC.")

    return parser.parse_args()

def validate_args(args):
    """
    ensure command line arguments are sane
    Args:
        args: (argparse.Namespace)
    Returns
        (str/None): if error a string describing problen else None
    """
    if not args.input.exists():
        return f"file {args.input} does not exist!"

    if args.resolution <= 0.0:
        return f"resolution must be greater than zero!"

    return None

def main():
    args = get_args()
    message = validate_args(args)

    if message is not None:
        print(f"Error {message}", file=sys.stderr, flush=True)
        return

    # ensure output has correct suffix
    out_file = args.output.with_suffix('.mrc')

    scale = refine_volume_data(args.input, out_file, args.resolution)

    print(f"{args.input} scaled to {scale} and written to {out_file}")

if __name__ == "__main__":
    main()
