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

Generating tetrahedral meshes from pixel data
"""
# set up linting
# pylint: disable = import-error

import sys
import argparse
import pathlib
from enum import Enum

from ffeamesh.sixtets import (convert_mrc_to_6tets, convert_mrc_to_6tets_interp)

class Method(Enum):
    """
    enumeration of possible methods of converting
    """
    ## simple decomp into 5 tets
    PLAIN = 'plain'

    ## interpolated decomp into 5 tets
    INTERP = 'interp'

def get_args():
    """
    get the command line arguments
    Returns:
        (argparse.namespace)
    """
    parser = argparse.ArgumentParser("""process MRC files and produces a regular
        tetrahedral volumetric mesh for FFEA using the "marching tet" algorithm.
        This is written out in the tetgen .ele, .face, and .node file format for
        later use in FFEA, and .vtk for mesh analysis.

        Coding:   Molly Gravett (bsmgr@leeds.ac.uk),
        Joanna Leng (J.Leng@leeds.ac.uk),
        Jarvellis Rogers (J.F.Rogers1@leeds.ac.uk)
        Jonathan Pickering (J.H.Pickering@leeds.ac.uk)""")

    parser.add_argument("-i",
                        "--input",
                        type=pathlib.Path,
                        required=True,
                        help="input file")

    parser.add_argument("-o",
                        "--output",
                        type=pathlib.Path,
                        required=True,
                        help="output file name root, (no suffix)")

    parser.add_argument("-v",
                        "--vtk",
                        action="store_true",
                        help="produce vtk output.")

    parser.add_argument("-f",
                        "--ffea",
                        action="store_true",
                        help="produce ffea output.")

    parser.add_argument("-t",
                        "--threshold",
                        type=float,
                        default=0.0,
                        help="""lower filter for voxels, default zero""")

    parser.add_argument("-w",
                        "--overwrite",
                        action="store_true",
                        help="overwrite")

    parser.add_argument("-m",
                        "--method",
                        type=Method,
                        default=Method.PLAIN,
                        help="choice of method used to make tetrohedrons",
                        choices=list(Method))

    parser.add_argument("-V",
                        "--verbose",
                        action="store_true",
                        help="write verbose output")

    parser.add_argument("-p",
                        "--progress",
                        action="store_true",
                        help="print progress during operation (may slow package)")

    return parser.parse_args()

def validate_command_line(args):
    """
    validate the command line arguments#
    Args:
        args (argeparse.Namespace): the command line arguments
    Returns
        (str): error message if fail else None
    """
    # check the input file can be found
    if not args.input.exists():
        return f"Error: file {args.input} does not exist!"

    # check for overwriteing files
    if not args.overwrite:
        if args.vtk:
            file = args.output.with_suffix(".vtk")
            if file.exists():
                return f"Error: file {file} exists, use option -w to allow overwrite."

        if args.ffea:
            files = []
            files.append(args.output.with_suffix(".vtk"))
            files.append(args.output.with_suffix(".1.ele"))
            files.append(args.output.with_suffix(".1.face"))
            files.append(args.output.with_suffix(".1.node"))
            for file in files:
                if file.exists():
                    return f"Error: file {file} exists, use option -w to allow overwrite."

    # check that some output has been specified
    if not args.vtk and not args.ffea:
        return "Error: you must specify and output type (vtk and/or ffea)"

    return None

def main():
    """
    run the script
    """
    args = get_args()
    error_message = validate_command_line(args)

    if error_message is not None:
        print(error_message, file=sys.stderr)
        return

    if args.method == Method.PLAIN:
        convert_mrc_to_6tets(args.input,
                             args.output,
                             args.threshold,
                             args.ffea,
                             args.vtk,
                             args.verbose,
                             args.progress)

    elif args.method == Method.INTERP:
        convert_mrc_to_6tets_interp(args.input,
                                    args.output,
                                    args.threshold,
                                    args.ffea,
                                    args.vtk,
                                    args.verbose,
                                    args.progress)


if __name__ == "__main__":
    main()
