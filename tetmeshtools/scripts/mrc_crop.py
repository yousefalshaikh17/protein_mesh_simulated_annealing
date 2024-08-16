#!/usr/bin/env python3
"""

 Script that crops 3D mrc image file data.

-----------------

Licensed under the GNU General Public License, Version 3.0 (the "License"); you
may not use this file except in compliance with the License. You may obtain a
copy of the License at <https://www.gnu.org/licenses/gpl-3.0.html>.

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

This work was funded by Joanna Leng's EPSRC funded RSE Fellowship (EP/R025819/1)

@copyright 2023
@author: j.h.pickering@leeds.ac.uk and j.leng@leeds.ac.uk
"""
# set up linting conditions
# pylint: disable = import-error

import argparse
import pathlib
import numpy as np
import mrcfile

def get_args():
    """
    get the command line arguments
        Returns
            (argparse.namespace)
    """
    parser = argparse.ArgumentParser(description="crops 3D mrc image file data")

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

    parser.add_argument("-y",
                        "--ystart",
                        type=int,
                        default=0,
                        help="start of x section")

    parser.add_argument("-Y",
                        "--yend",
                        type=int,
                        default=100,
                        help="end of x section")

    parser.add_argument("-z",
                        "--zstart",
                        type=int,
                        default=0,
                        help="start of x section")

    parser.add_argument("-Z",
                        "--zend",
                        type=int,
                        default=100,
                        help="end of x section")

    parser.add_argument("-p",
                        "--pad",
                        action='store_true',
                        help="add zero padding one voxel thick around selection")

    return parser.parse_args()

def output_mrcfile(data, out_file, original_mrc):
    """
    write out a mrc file
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
    """
    slice out the required part of the data file
    Args:
        mrc (mrcfile.file): source file
        args (argpars.namespace): the command line arguments
    """
    data = mrc.data[args.zstart:args.zend, args.ystart:args.yend, args.xstart:args.xend]
    print(f"{type(data)}, {data.shape}")
    if args.pad:
        data = np.pad(data, 1)
    output_mrcfile(data, args.output, mrc)

def main():
    """
    run the script
    """
    args = get_args()
    with mrcfile.open(args.input, mode='r+') as mrc:
        chop_image(mrc, args)

    print(f"{args.input} sliced and written to {args.output}")

if __name__ == "__main__":
    main()
