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

 Created on Mon Oct 25 11:50:39 2021

 @author: mollygravett
 modified jonathan pickering 23Aug22
"""
# set up linting
# pylint: disable = import-error

import mrcfile
import numpy as np
import ffeamesh.mrclininterpolate as mi

def scale_int_list(int_list, scale_factor):
    """
    scale an int list by a float factor convert the results
    to ints and ensure all are greater than zero
    Args:
        int_list ([int]): the list to be scaled
        scale_factor (float): the factor
    Returns:
        ([int]): scaled list
    Raises:
        ValueError is any entry is zero
    """
    out_list = [int(x*scale_factor) for x in int_list]
    if not all(x>0.0 for x in out_list):
        raise ValueError("Scale factor too small, one axis has no voxels!")
    return out_list

def make_new_density_array(vox_counts, header, image):
    """
    make a numpy density array that can be used in a mrc file
    Args:
        voxel_counts [int]: the number of voxes on x, y & z axis
    Retruns
        (numpy.array of float)
    """
    array = np.zeros(vox_counts, dtype=np.float32)

    step_x = header.cella.x/array.shape[2]
    step_y = header.cella.y/array.shape[1]
    step_z = header.cella.z/array.shape[0]

    start_x = header.origin.x + (step_x/2.0)
    start_y = header.origin.y + (step_y/2.0)
    start_z = header.origin.z + (step_z/2.0)

    # iterate the voxels and make tet connectivities
    for voxel_z in range(vox_counts[0]):
        for voxel_y in range(vox_counts[1]):
            for voxel_x in range(vox_counts[2]):
                # if voxel_x > 10:
                #     array[voxel_z, voxel_y, voxel_x] = 1.0
                ctr_x = start_x + (voxel_x*step_x)
                ctr_y = start_y + (voxel_y*step_y)
                ctr_z = start_z + (voxel_z*step_z)
                density, distance = image.density_or_distance_at(ctr_x, ctr_y, ctr_z)
                if distance > 0.0:
                    raise ValueError("Attempt to sample density out of image")
                array[voxel_z, voxel_y, voxel_x] = density

    return array

def resample_volume_data(in_file, out_file, scale_factor):
    """
    resample the input file producing an ouput with a different number of voxels
    Args:
        in_file (pathlib.Path) the input file
        out_file (pathlib.Path) the output file
        scale_factor (float): range [0, 1], multipler for the number of voxels on each axis
    """
    with mrcfile.mmap(in_file, mode='r+') as mrc:
        image = mi.MRCImage(mrc)
        vox_counts = scale_int_list([image.get_nz(), image.get_ny(), image.get_nx()], scale_factor)
        array = make_new_density_array(vox_counts, mrc.header, image)
        write_array_to_mrcfil(array, out_file, mrc.header)

def write_array_to_mrcfil(array, out_file, header):
    """
    write out a scaled mrc file
        Args:
            array (np.array): array of voxels
            out_file (pathlib.Path) the output file
            header (mrc.header) header from original mrc file
    """
    print(f"Origin {header.origin}")
    print(f"array shape {array.shape}")

    with mrcfile.new(out_file, overwrite=True) as mrc:
        mrc.set_data(array)
        mrc.header.origin = header.origin
        mrc.header.cella = header.cella
        mrc.header.cellb = header.cellb
