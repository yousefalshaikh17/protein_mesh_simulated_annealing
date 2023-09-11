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

    Authors: Joanna Leng, Jonathan Pickering - University of Leeds
    Emails: J.Leng@leeds.ac.uk, J.H.Pickering@leeds.ac.uk
"""
import argparse
import pathlib
import mrcfile
import numpy as np

import ffeamesh.mrclininterpolate as mi

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
    description  = ("run a density scan across an mrc file on "
                    "x, y and z axis printing out locations and "
                    "image densities")
    parser = argparse.ArgumentParser(description)

    parser.add_argument("-i",
                        "--image_file",
                        type=mrc_file,
                        required=True,
                        help="mrc image file")

    return parser.parse_args()

def test_error(image_file):
    """
    test if the interpolation throws an exception
    Args:
        image_file: MRCImage
    Returns:
        number of exceptions raised
    """
    with mrcfile.open(image_file, mode='r+') as mrc:
        image = mi.MRCImage(mrc)
        values = np.linspace(-2.5, 27.5, num=100)
        #values = [16.0, 19.0, 21.0, 23.0]
        exceptions = 0

        for coord in values:
            try:
                image.density_or_distance_at(coord, 12.5, 12.5)
                image.density_or_distance_at(12.5, coord, 12.5)
                image.density_or_distance_at(12.5, 12.5, coord)
            except Exception as exc:
                print(f">>> failed {coord} {exc}")
                exceptions += 1

    print(f"Test run {exceptions} exceptions thrown")

    return exceptions

def print_interpolation(image_file):
    """
    print results of interpolation
    Args:
        image_file: MRCImage
    """
    with mrcfile.open(image_file, mode='r+') as mrc:
        image = mi.MRCImage(mrc)
        values = np.linspace(-2.5, 27.5, num=100)
        #values = [16.0, 19.0, 21.0, 23.0]
        # file = pathlib.Path("myout.csv")
        # of = file.open('w')
        for coord in values:
            found_density, dist2 = image.density_or_distance_at(coord, 12.5, 12.5)
            # print(f"{coord}, {found_density}, {dist2}", file=of)
            print(f"({coord}, 12.5, 12.5) => {found_density}, {dist2}")
            found_density, dist2 = image.density_or_distance_at(12.5, coord, 12.5)
            print(f"(12.5, {coord}, 12.5) => {found_density}, {dist2}")
            found_density, dist2 = image.density_or_distance_at(12.5, 12.5, coord)
            print(f"(12.5, 12.5, {coord}) => {found_density}, {dist2}\n")
        # of.close()

def main():
    """
    tests the interplation function across an image file
    """
    args = get_args()
    if test_error(args.image_file) == 0:
        print_interpolation(args.image_file)

if __name__ == "__main__":
    main()
