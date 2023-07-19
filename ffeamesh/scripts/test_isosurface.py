#!/usr/bin/env python3
"""
 mrc_crop.pylint

 Script that crops 3D mrc image file data.

-----------------

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
import numpy as np
import mrcfile

import ffeamesh.tetmeshtools.tetgenread as tr
import ffeamesh.optimizemesh.costfunction as cf
import ffeamesh.mrclininterpolate as mi

def parse_position(text):
    """
    convert text to position
    Args:
        text (str): user input
    Returns
        [float]
    """
    position = [float(part) for part in text.split()]

    if len(position) != 3:
        raise ValueError()

    return position

def tetgen_root_file(name):
    """
    check string is tetgen root file name
    Args:
        name str
    Returns:
        pathlib.Path
    """
    if not name.endswith(".1"):
        raise ValueError(f"File {name} is not tetgen root file")

    for suffix in [".node", ".ele", ".face"]:
        file = pathlib.Path(name + suffix)
        if not file.exists():
            raise ValueError(f"File {file} does not exist!")

    return pathlib.Path(name)

def mrc_file(name):
    """
    check string is mrc image file
    Args:
        name (str)
    Returns
        pathlib.Path
    """
    path = pathlib.Path(name)

    if not path.suffix == ".mrc":
        raise ValueError()

    if not path.exists():
        raise ValueError()

    return path

def add_file_arguments(s_file):
    """
    add args for 'file' option
    Args:
        s_file (subparser)
    """
    s_file.add_argument("-f",
                        "--input_file",
                        type=tetgen_root_file,
                        required=True,
                        help="input file")

    s_file.add_argument("-i",
                        "--image_file",
                        type=mrc_file,
                        required=True,
                        help="mrc image file")

    s_file.add_argument("-o",
                        "--output",
                        type=pathlib.Path,
                        required=True,
                        help="output file")

    s_file.add_argument("-v",
                        "--isovalue",
                        required=True,
                        type=float,
                        help="the isovalue")

def add_pos_arguments(sp_pos):
    """
    add args for 'file' option
    Args:
        s_file (subparser)
    """
    sp_pos.add_argument("-i",
                        "--image_file",
                        type=mrc_file,
                        required=True,
                        help="mrc image file")

    sp_pos.add_argument("-p",
                        "--position",
                        required=True,
                        type=float,
                        nargs='+',
                        help="the position")

def get_args():
    """
    get the command line arguments
        Returns
            (argparse.namespace)
    """
    parser = argparse.ArgumentParser("""test isosurface""")
    subparser = parser.add_subparsers(dest="subparser")

    s_file = subparser.add_parser("file",
                                  help="read tetgen format file and run")

    sp_pos = subparser.add_parser("position",
                                  help="calculat at one position")
    add_file_arguments(s_file)
    add_pos_arguments(sp_pos)

    return parser.parse_args()

def main():
    """run the program"""
    args = get_args()

    if args.subparser == "file":
        run_file(args)
    elif args.subparser == "position":
        run_position(args)
    return

def run_file(args):
    """
    run on a tetgen file
    Args:
        args (argparse.namespace)
    """
    model = tr.make_model_from_tetgen(args.input_file)
    with mrcfile.open(args.image_file, mode='r+') as mrc:
        image = mi.MRCImage(mrc)
        total_density, total_distance = cf.fit_to_isovalue(model.get_surface(),
                                                           args.isovalue,
                                                           image)

        print(total_density, total_distance)

def run_position(args):
    """
    calculat for a position
    Args:
        args (argparse.namespace)
    """
    with mrcfile.open(args.image_file, mode='r+') as mrc:
        image = mi.MRCImage(mrc)
        image.depress_below_isovalue(0.1)
        found_density, dist2 = image.density_or_distance_at(args.position[0],
                                                            args.position[1],
                                                            args.position[2])

        print(found_density, dist2)

if __name__ == "__main__":
    main()
