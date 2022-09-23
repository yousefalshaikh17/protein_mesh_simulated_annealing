# -*- coding: utf-8 -*-
#
#  This file is part of the FFEA simulation package
#
#  Copyright (c) by the Theory and Development FFEA teams,
#  as they appear in the README.md file.
#
#  FFEA is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  FFEA is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
#
#  To help us fund FFEA development, we humbly ask that you cite
#  the research papers on the package.
#

"""
Authors: Joanna Leng, Jonathan Pickering - University of Leeds
Emails: J.Leng@leeds.ac.uk, J.H.Pickering@leeds.ac.uk
"""
import sys
import pathlib
import argparse
from enum import Enum
from datetime import datetime
import getpass
import numpy as np

from ffeamesh.mrc_utility import (CellSize, CellAngles, CellProps, write_mrcfile)

class Voxtest(Enum):
    """
    enumeration of the three built in tests
    """
    TEST0 = "test0"
    TEST1 = "test1"
    TEST2 = "test2"
    TEST3 = "test3"

def get_args():
    """
    get the command line arguments
        Returns
            (argparse.namespace)
    """
    parser = argparse.ArgumentParser(description="""make simple mrc image files
        for use in testing""")

    parser.add_argument("-o",
                        "--output",
                        type=pathlib.Path,
                        required=True,
                        help="output file")

    parser.add_argument("-v",
                        "--voxtest",
                        type=Voxtest,
                        required=True,
                        help="which test do you require",
                        choices=list(Voxtest))

    parser.add_argument("-w",
                        "--overwrite",
                        action="store_true",
                        help="if output file exists overwrite")

    return parser.parse_args()

def validate_args(args):
    """
    check the arguments are applicable
        Args:
            args (argparse.Namespace): the command line arguments
        Returns:
            (str): error message, or None if no problem
    """
    if args.output.exists() and not args.overwrite:
        return f"Error input file {args.output} does exists, use -w to force overwrite."

    return None

def voxel_test(voxtest):
    """
    make data for built in tests
        Args:
            (Voxtest): enum specifing the selected test
        Returns:
            (numpy.ndarray data=np.float32), (CellProps) data plus cell size & angles
    """
    mrc_data = None

    if voxtest == Voxtest.TEST0:
        vox_data = [[[1,2],[3,4]],[[2,3],[4,5]]] # List nested as z > y > x
        mrc_data = np.tile(vox_data,1).astype(np.float32)

    elif voxtest == Voxtest.TEST1:
        vox_data = [[[2,3],[1,2]],[[1,2],[1,1]],[[2,3],[1,2]]]
        mrc_data = np.tile(vox_data,1).astype(np.float32)

    elif voxtest == Voxtest.TEST2:
        vox_data = [[[3,3],[2,2],[2,2],[3,3]],
                    [[2,2],[1,1],[1,1],[2,2]],
                    [[2,2],[1,1],[1,1],[2,2]],
                    [[3,3],[2,2],[2,2],[3,3]]]
        mrc_data = np.tile(vox_data,1).astype(np.float32)

    elif voxtest == Voxtest.TEST3:
        vox_data = [[[1]]]
        mrc_data = np.tile(vox_data,1).astype(np.float32)

    props = CellProps(CellSize(20.0, 20.0, 20.0), CellAngles(90.0, 90.0, 90.0))

    return mrc_data, props

def main():
    """
    run the script
    """
    args = get_args()
    error_message = validate_args(args)

    if error_message is not None:
        print(error_message, file=sys.stderr)
        return

    if args.voxtest is not None:
        test_image, cell_size = voxel_test(args.voxtest)
        label = str(args.voxtest)
        label += f" by {getpass.getuser()} on "
        label += datetime.now().strftime("%d-%b-%Y (%H:%M:%S.%f)")
        write_mrcfile(test_image, cell_size, args.output, label, args.overwrite)
        return

if __name__ == "__main__":
    main()
