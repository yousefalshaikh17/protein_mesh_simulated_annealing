#!/usr/bin/env python
"""
 A script that produces lower resolution mrc files.

 ----------------------------

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
