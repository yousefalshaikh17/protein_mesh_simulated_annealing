#!/usr/bin/env python
"""
 mrc_coarsen.py

 A script that produces lower resolution mrc files.

 ----------------------------

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

"""
# set up linting
# pylint: disable = import-error
import sys
import argparse
import pathlib
import mrcfile
import numpy as np

import tetmeshtools.mrclininterpolate as mi

from tetmeshtools.mrc_utility import voxel_size

def check_counts(in_str):
    """
    check the numbers of voxels are greter than zero
    """
    my_int = int(in_str)
    if my_int < 1:
        raise ValueError("number of voxels must be greater than zero")
    return my_int

def get_args():
    """
    run argparse to get command line
    Returns
        namspace
    """
    parser = argparse.ArgumentParser(description=" produces lower resolution mrc files")


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

    parser.add_argument("-n",
                        "--voxelcounts",
                        type=check_counts,
                        nargs='+',
                        required=True,
                        help="number of voxels on x y z axis")

    return parser.parse_args()

def write_array_to_mrcfile(array, out_file, header, origin, cell_size, voxel_counts):
    """
    write out a scaled mrc file
        Args:
            array (np.array): array of voxels
            out_file (pathlib.Path) the output file
            header (mrc.header) header from original mrc file
    """
    with mrcfile.new(out_file, overwrite=True) as mrc:
        mrc.set_data(array)
        mrc.header.origin.x = origin[0]
        mrc.header.origin.y = origin[1]
        mrc.header.origin.z = origin[2]

        mrc.header.cella.x = cell_size[0]
        mrc.header.cella.y = cell_size[1]
        mrc.header.cella.z = cell_size[2]

        mrc.header.nx = voxel_counts[0]
        mrc.header.mx = voxel_counts[0]

        mrc.header.ny = voxel_counts[1]
        mrc.header.my = voxel_counts[1]

        mrc.header.nz = voxel_counts[2]
        mrc.header.mz = voxel_counts[2]

        mrc.header.cellb = header.cellb

def make_data_array(image, voxel_counts, new_vox_size, new_origin):
    """
    make the main data array
    Args:
        image (MRCImage)
        voxel_counts [int int int]: number of voxels on x, y & z axis
        new_vox_size [float, float, float]: size of voxel x, y, z
        new_origin [float, float, float]: location of cell
    Returns
        numpy.ndarray the densities
    """
    half_new_vox_size = [x/2.0 for x in new_vox_size]

    #array = np.zeros(voxel_counts, dtype=np.float32)
    array = np.zeros([voxel_counts[2], voxel_counts[1], voxel_counts[0]], dtype=np.float32)
    for z_index in range(0, voxel_counts[2]):
        z_coord = (float(z_index)*new_vox_size[2]) + new_origin[2] + half_new_vox_size[2]

        for y_index in range(0, voxel_counts[1]):
            y_coord = (float(y_index)*new_vox_size[1]) + new_origin[1] + half_new_vox_size[1]

            for x_index in range(0, voxel_counts[0]):
                x_coord = (float(x_index)*new_vox_size[0]) + new_origin[0] + half_new_vox_size[0]

                density, distance = image.density_or_distance_at(x_coord, y_coord, z_coord)
                if distance != 0.0:
                    coord = f"({x_coord}, {y_coord}, {z_coord})"
                    index = f"({x_index}, {y_index}, {z_index})"
                    raise ValueError(f"Point {coord} {index} outside image {distance}")

                array[z_index, y_index, x_index] = density
                #print(f"({x_coord:7.2f}, {y_coord:7.2f}, {z_coord:7.2f}) => {density}")

    return array

def main(input_file, output_file, voxel_counts):
    """
    do the conversion
    """
    with mrcfile.mmap(input_file, mode='r+') as mrc:
        image = mi.MRCImage(mrc)

        # get vooxel sizes
        delta = voxel_size(mrc)

        # find new coordinate origin
        new_origin = []
        new_origin.append(np.float32(mrc.header.origin['x'])+(delta.dx/2.0))
        new_origin.append(np.float32(mrc.header.origin['y'])+(delta.dy/2.0))
        new_origin.append(np.float32(mrc.header.origin['z'])+(delta.dz/2.0))

        # new cell size
        new_cell_size = []
        new_cell_size.append(mrc.header.cella.x-delta.dx)
        new_cell_size.append(mrc.header.cella.y-delta.dy)
        new_cell_size.append(mrc.header.cella.z-delta.dz)

        # new voxel size
        new_vox_size = []
        new_vox_size.append(new_cell_size[0]/voxel_counts[0])
        new_vox_size.append(new_cell_size[1]/voxel_counts[1])
        new_vox_size.append(new_cell_size[2]/voxel_counts[2])

        print(voxel_counts)
        print(new_origin)
        print(new_cell_size)
        print(new_vox_size)

        array = make_data_array(image,
                                voxel_counts,
                                new_vox_size,
                                new_origin)

        write_array_to_mrcfile(array,
                               output_file,
                               mrc.header,
                               new_origin,
                               new_cell_size,
                               voxel_counts)

if __name__ == "__main__":
    args = get_args()
    if len(args.voxelcounts) != 3:
        sys.exit("You must provide exactly 3 voxel counts")

    main(args.input, args.output, args.voxelcounts)
