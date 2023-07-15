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

def get_args():
    """
    get the command line arguments
        Returns
            (argparse.namespace)
    """
    parser = argparse.ArgumentParser("""test isosurface""")

    parser.add_argument("-f",
                        "--input_file",
                        type=tetgen_root_file,
                        required=True,
                        help="input file")

    parser.add_argument("-i",
                        "--image_file",
                        type=mrc_file,
                        required=True,
                        help="mrc image file")

    parser.add_argument("-o",
                        "--output",
                        type=pathlib.Path,
                        required=True,
                        help="output file")

    parser.add_argument("-v",
                        "--isovalue",
                        required=True,
                        type=float,
                        help="the isovalue")

    return parser.parse_args()

def main():
    """run the program"""
    args = get_args()

    model = tr.make_model_from_tetgen(args.input_file)
    with mrcfile.open(args.image_file, mode='r+') as mrc:
        image = mi.MRCImage(mrc)
        total_density, total_distance = cf.fit_to_isovalue(model.get_surface(),
                                                           args.isovalue,
                                                           image)

        print(total_density, total_distance)

if __name__ == "__main__":
    main()
