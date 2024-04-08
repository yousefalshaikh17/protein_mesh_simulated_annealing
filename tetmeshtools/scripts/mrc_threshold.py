#!/usr/bin/env python3
"""
  A script that reads a mrc file and set all voxels less than than the
  thershold to zero, then outputs to a new file.

  ----------------------

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
import sys
import datetime
import getpass
import mrcfile
from tetmeshtools.mrc_utility import (threshold_mrc_image,
                                  get_cell_props,
                                  write_mrcfile)

def get_args():
    """
    get the command line arguments
        Returns
            (argparse.namespace)
    """
    description = ("reads a mrc file and set all voxels less than than "
                   "the thershold to zero, then outputs to a new file")
    parser = argparse.ArgumentParser(description=description)

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

    parser.add_argument("-t",
                        "--threshold",
                        type=float,
                        required=True,
                        help="""filter limit all voxels with values less
                                than or equal will be set to zero""")

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

    if not args.overwrite and args.output.exists():
        return f"file {args.output} already exists run with -w to overwrite"

    return None

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
        data = threshold_mrc_image(mrc.data, args.threshold)
        date = datetime.datetime.now().strftime("%x")
        label = f"made by {getpass.getuser()}, on {date}"
        write_mrcfile(data, get_cell_props(mrc), args.output, label, args.overwrite)
        print(f"thresholded output from {args.input} written to {args.output}")

if __name__ == "__main__":
    main()
