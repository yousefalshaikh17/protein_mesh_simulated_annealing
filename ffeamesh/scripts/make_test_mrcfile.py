#!/usr/bin/env python3
"""
 make_test+mrcfile.pylint

 Script that makes simple mrc image files for use in testing.

 -------------------------------

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

import sys
import pathlib
import argparse
from enum import Enum
from datetime import datetime
import getpass
import numpy as np

from ffeamesh.mrc_utility import (CellSize, CellAngles, CellProps, write_mrcfile)
from ffeamesh.comlinesupport import (positive_int, positive_float)

class Voxtest(Enum):
    """
    enumeration of the three built in tests
    """
    CUBE = "cube"
    SPHERE = "sphere"
    DBELL = "dbell"
    SIMPLE = "simple"

def get_args():
    """
    get the command line arguments
        Returns
            (argparse.namespace)
    """
    description = "makes simple mrc image file for use in testing"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("-o",
                        "--output",
                        type=pathlib.Path,
                        required=True,
                        help="output file")

    parser.add_argument("-v",
                        "--voxtest",
                        type=Voxtest,
                        required=True,
                        help=f"which test do you require, {[el.value for el in Voxtest]}")

    parser.add_argument("-c",
                        "--voxel_count",
                        type=positive_int,
                        default=10,
                        help="the number of voxels per axis (total = n*n*n)")

    parser.add_argument("-l",
                        "--voxel_size",
                        type=positive_float,
                        default=5.0,
                        help="size of voxel")

    parser.add_argument("-s",
                        "--object_size",
                        type=positive_float,
                        default=5.0,
                        help="size of object")

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

def voxel_test(voxtest, count, voxel_size, object_size):
    """
    make data for built in tests
        Args:
            (Voxtest): enum specifing the selected test
            count (int): number of cells on each dimension
            voxel_size (float): size of voxel
            object_size (float): size of object
        Returns:
            (numpy.ndarray data=np.float32), (CellProps) data plus cell size & angles
    """
    mrc_data = None
    cell_size = voxel_size * count

    if voxtest == Voxtest.CUBE:
        mrc_data = make_cube(count)

    elif voxtest == Voxtest.SPHERE:
        mrc_data = make_sphere(voxel_size, object_size, count)

    elif voxtest == Voxtest.DBELL:
        mrc_data = make_dbell(voxel_size, object_size, count)

    elif voxtest == Voxtest.SIMPLE:
        mrc_data = make_simple()

    props = CellProps(CellSize(cell_size, cell_size, cell_size),
                      CellAngles(90.0, 90.0, 90.0))

    return mrc_data, props

def make_simple():
    """
    make a simple gradient
    Returns
        np array
    """
    test_image = np.full((5, 5, 5), dtype=np.float32, fill_value=0.2)
    fill_test(test_image)

    return test_image

def fill_test(image):
    """
    fill the image array
    Args:
        image: numpy.ndarray
    """
    for i in range(1, 4):
        for j in range(1, 4):
            for k in range(1, 4):
                image[i][j][k] = 0.8

    image[2][2][2] = 1.2

def make_cube(count):
    """
    make the array for a cube
    Args:
        voxel_size (float): voxel dimension
        count (int): number of cells on each dimension
    Returns
        np array
    """
    vox_data = np.full((count, count, count), 0.05, dtype=np.float32)

    half_count = round(count/2)

    fill_cube(vox_data, half_count, 5, 0.1)
    fill_cube(vox_data, half_count, 4, 0.2)
    fill_cube(vox_data, half_count, 3, 0.4)
    fill_cube(vox_data, half_count, 2, 0.8)
    fill_cube(vox_data, half_count, 1, 1.6)

    vox_data[half_count, half_count, half_count] = 1.0

    return vox_data

def fill_cube(vox_data, half_count, offset, level):
    """
    fill a cube of data
    vox_data (numpy array): the array
    half_count (int): half the array size
    offset (int): the range around the ctr to be set to level
    level (float): the value to be set
    """
    for i in range(half_count-offset, half_count+offset):
        for j in range(half_count-offset, half_count+offset):
            for k in range(half_count-offset, half_count+offset):
                vox_data[i, j, k] = level

def make_dbell(voxel_size, object_size, count):
    """
    make the array for a sphere
    Args:
        voxel_size (float): voxel dimension
        count (int): number of cells on each dimension
    Returns
        np array
    """
    vox_data = np.full((count, count, count), 0.05, dtype=np.float32)

    frac_centre = (count * voxel_size)/10.0
    centre = frac_centre*5

    fill_sphere(count,
                (centre, centre, 3*frac_centre),
                object_size,
                voxel_size,
                vox_data)
    fill_sphere(count,
                (centre, centre, 7*frac_centre),
                object_size,
                voxel_size,
                vox_data)

    return vox_data

def make_sphere(voxel_size, object_size, count):
    """
    make the array for a sphere
    Args:
        voxel_size (float): cell dimension
        count (int): number of cells on each dimension
        object_size (float): size of sphere
    Returns
        np array
    """
    vox_data = np.full((count, count, count), 0.05, dtype=np.float32)

    centre = (count * voxel_size)/2.0

    fill_sphere(count,
                (centre, centre, centre),
                object_size,
                voxel_size,
                vox_data)

    return vox_data

def fill_sphere(count, centre, object_size, voxel_size, vox_data):
    """
    fill an array with spherical data densities
    Args:
        count (int): number of cells on each dimension
        centre ((float)): cartesian coordinates of sphere centre
        object_size (float): the size of the sphere
        voxel_size (float): the size of a voxel
        vox_data (np.array(float)): the image array
    Returns
        np.array
    """
    half_size = voxel_size/2.0
    for i in range(count):
        del_z = np.float32(i)*voxel_size + half_size
        del_z = centre[2] - del_z
        for j in range(count):
            del_y = np.float32(j)*voxel_size + half_size
            del_y = centre[1] - del_y
            for k in range(count):
                del_x = np.float32(k)*voxel_size + half_size
                del_x = centre[0] - del_x
                radius = np.sqrt(del_x*del_x + del_y*del_y + del_z*del_z)

                if radius < object_size*0.5:
                    vox_data[i, j, k] += 1.0
                elif radius < object_size:
                    vox_data[i, j, k] += 0.75
                elif radius < object_size*1.5:
                    vox_data[i, j, k] += 0.5
                elif radius < object_size*2:
                    vox_data[i, j, k] += 0.25
                elif radius < object_size*2.5:
                    vox_data[i, j, k] += 0.125

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
        test_image, voxel_size = voxel_test(args.voxtest,
                                            args.voxel_count,
                                            args.voxel_size,
                                            args.object_size)
        label = str(args.voxtest)
        label += f" by {getpass.getuser()} on "
        label += datetime.now().strftime("%d-%b-%Y (%H:%M:%S.%f)")
        write_mrcfile(test_image, voxel_size, args.output, label, args.overwrite)
        return

if __name__ == "__main__":
    main()
