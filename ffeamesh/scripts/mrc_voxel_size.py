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

@author: jonathan pickering 23Aug22
"""
# set up linting conditions
# pylint: disable = import-error

import argparse
import pathlib
import mrcfile

from ffeamesh.mrc_utility import voxel_size

def get_args():
    """
    get the command line arguments
        Returns
            (argparse.namespace)
    """
    parser = argparse.ArgumentParser(description="print an MRC file's voxel size")

    parser.add_argument("-i",
                        "--input",
                        type=pathlib.Path,
                        required=True,
                        help="input file")

    return parser.parse_args()

def main():
    """
    run the script
    """
    args = get_args()
    with mrcfile.open(args.input, mode='r+') as mrc:
        delta = voxel_size(mrc)

    print(f"{args.input} voxel size {round(delta.dx, 2)}, "\
          f"{round(delta.dy, 2)}, {round(delta.dz, 2)}")

if __name__ == "__main__":
    main()
