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
import argparse
import pathlib

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

def main():
    args = get_args()

    if not args.input.exists():
        print(f"Input file {args.input} does not exist!")
        return

    out_file = args.output
    if not out_file.suffix == '.mrc':
        out_file = out_file.with_suffix('.mrc')

    scale = refine_volume_data(args.input, out_file, args.resolution)

    print(f"{args.input} scaled to {scale} and written to {out_file}")

if __name__ == "__main__":
    main()
