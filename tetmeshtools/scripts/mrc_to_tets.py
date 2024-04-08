#!/usr/bin/env python
"""
 five_tets.py

 A script that processes MRC files and produces a regular
 tetrahedral volumetric mesh, with five or six tetrahedra
 in each voxel, using the "marching tet" algorithm. This is
 written out in the tetgen .ele, .face, and .node file format,
 and .vtk for mesh analysis.

 ----------------------------

Licensed under the GNU General Public License, Version 3.0 (the "License"); you
may not use this file except in compliance with the License. You may obtain a
copy of the License at <https://www.gnu.org/licenses/gpl-3.0.html>.

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

This work was funded by Joanna Leng's EPSRC funded RSE Fellowship (EP/R025819/1)

@copyright 2020
@author: j.h.pickering@leeds.ac.uk and j.leng@leeds.ac.uk

"""
# set up linting
# pylint: disable = import-error

import sys
import argparse
import pathlib
from tetmeshtools.maketets import convert_mrc_to_tets

def get_args():
    """
    get the command line arguments
    Returns:
        (argparse.namespace)
    """
    description = ("process MRC files and produces a regular "
        "tetrahedral volumetric mesh using the \"marching tet\" algorithm. "
        "Output in the tetgen .ele, .face, and .node file format, and .vtk for mesh analysis."

        "Coding:   Molly Gravett (bsmgr@leeds.ac.uk), "
        "Joanna Leng (J.Leng@leeds.ac.uk), "
        "Jarvellis Rogers (J.F.Rogers1@leeds.ac.uk), "
        "Jonathan Pickering (J.H.Pickering@leeds.ac.uk)")

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
                        help="output file name root, (no suffix)")

    parser.add_argument("-v",
                        "--vtk",
                        action="store_true",
                        help="output files in vtk format.")

    parser.add_argument("-f",
                        "--ftetg",
                        action="store_true",
                        help="output files in tetgen format.")

    parser.add_argument("-t",
                        "--threshold",
                        type=float,
                        default=0.0,
                        help="""lower filter for voxels, default zero""")

    parser.add_argument("-w",
                        "--overwrite",
                        action="store_true",
                        help="overwrite")

    parser.add_argument("-V",
                        "--verbose",
                        action="store_true",
                        help="write verbose output")

    parser.add_argument("-p",
                        "--progress",
                        action="store_true",
                        help="print progress during operation (may slow package)")

    parser.add_argument("-6",
                        "--use_six_tets",
                        action="store_true",
                        help="decompose into six tets, default 5")

    parser.add_argument("-m",
                        '--low_vertices',
                        type=int,
                        default=2,
                        choices=[1, 2, 3, 4],
                        help="number vertices < isovalue for tet cull (default: %(default)s)")

    parser.add_argument("-n",
                        "--vox_counts",
                        type=int,
                        nargs='+',
                        help="number of voxels on x y z axis")

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

        if args.ftetg:
            files = []
            files.append(args.output.with_suffix(".vtk"))
            files.append(args.output.with_suffix(".1.ele"))
            files.append(args.output.with_suffix(".1.face"))
            files.append(args.output.with_suffix(".1.node"))
            for file in files:
                if file.exists():
                    return f"Error: file {file} exists, use option -w to allow overwrite."

    # check that some output has been specified
    if not args.vtk and not args.ftetg:
        return "Error: you must specify and output type (vtk and/or ftetg)"

    if args.vox_counts is not None:
        if len(args.vox_counts) != 3:
            return f"Voxel counts must have 3 components, {args.vox_counts}"

        if not all(x>1 for x in args.vox_counts):
            return f"Voxel counts must be > 1, {args.vox_counts}"

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

    # pr = cProfile.Profile()
    # pr.enable()
    convert_mrc_to_tets(args.input,
                        args.output,
                        args.threshold,
                        args.ftetg,
                        args.vtk,
                        args.verbose,
                        args.progress,
                        args.vox_counts,
                        args.low_vertices,
                        args.use_six_tets)
    # pr.disable()
    # s = io.StringIO()
    # sortby = pstats.SortKey.CUMULATIVE
    # ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    # ps.print_stats(10)
    # print(s.getvalue())

if __name__ == "__main__":
    main()
