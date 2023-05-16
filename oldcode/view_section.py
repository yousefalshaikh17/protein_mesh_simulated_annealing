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
# set up linting
# pylint: disable = import-error

import argparse
import pathlib
import itertools
import numpy as np
import mrcfile
from PIL import Image

import ffeamesh.mrclininterpolate as mli

def display_cross_section(image, xsize, zsize, y_coord, isovalue, file_name):
    """
    display a section of the file
    Args:
        image (mrcfile): image to analyse
        xsize (float): plane size in x
        zsize (float): plane size in z
        y_coord (float): y coordinate of slice
        isovalue (float): the target isovalue
        file_name (str): the root name of the output file
    """
    x_range = np.linspace(image.low_limit_x,
                          image.high_limit_x-image.offset_x*2,
                          xsize)
    z_range = np.linspace(image.low_limit_z,
                          (image.high_limit_z-image.offset_z*2),
                          zsize)

    y_coord += image.low_limit_y

    image_slice = np.zeros((len(x_range), len(z_range)))

    for i, x_coord in enumerate(x_range):
        for j, z_coord in enumerate(z_range):
            image_slice[i, j], _ = image.density_or_distance_at(x_coord, y_coord, z_coord)

    show_image(image_slice, isovalue, file_name)

def show_image(slice_array, isovalue, file_name):
    """
    display the image slice as image and save it to file
    Args:
        slice_array (numpy.array(float)): 2d slice of 3d image
        isovalue (float): the target isovalue
        file_name (str): the root name of the output file
    """
    image_array = np.zeros(slice_array.shape, dtype=np.uint8)

    for i in range(image_array.shape[0]):
        for j in range(image_array.shape[1]):
            if slice_array[i, j] > isovalue:
                image_array[i, j] = 255

    img = Image.fromarray(image_array, "L")
    # Display the Numpy array as Image
    #img.show()
    img.save(file_name)

def get_args():
    """
    get command line
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("-i",
                        "--input",
                        type=pathlib.Path,
                        required=True,
                        help="MRC source file")

    parser.add_argument("-x",
                        "--xsize",
                        type=int,
                        required=True,
                        help="x steps")

    parser.add_argument("-z",
                        "--zsize",
                        type=int,
                        required=True,
                        help="z steps")

    parser.add_argument("-v",
                        "--isovalue",
                        type=float,
                        required=True,
                        help="value defining surface")

    return parser.parse_args()

def main():
    """
    demo of linear interpolation in data from mrc file:
    map_equilibrium_mystructrigor_15A_0p00202.mrc
    """
    args = get_args()
    #with mrcfile.open(args.input, mode='r+') as mrc:
    with mrcfile.mrcmemmap.MrcMemmap(args.input) as mrc:
        image = mli.MRCImage(mrc)
        count = itertools.count(1)
        total = 305-45
        for i in range(45, 305, 1):
            y_coord = float(i)
            file_name = f"yscan_{i}.png"
            display_cross_section(image, args.xsize, args.zsize, y_coord, args.isovalue, file_name)
            print(f"Image {next(count)} of {total}")

if __name__ == "__main__":
    main()
