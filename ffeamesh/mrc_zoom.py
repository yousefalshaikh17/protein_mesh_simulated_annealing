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
from scipy import ndimage
import ffeamesh.mrclininterpolate as mi
from ffeamesh.maketets import make_density_at_vertex_grid

def resample_volume_data(in_file, out_file, resolution):
    """
    resample the input file producing an ouput with a different number of voxels
    Args:
        in_file (pathlib.Path) the input file
        out_file (pathlib.Path) the output file
        resolution (float) the new resolution
    Return
    """
    with mrcfile.mmap(in_file, mode='r+') as mrc:
        image = mi.MRCImage(mrc)
        vox_counts = [image.get_nx(), image.get_ny(), image.get_nz()]//resolution
        grid = make_density_at_vertex_grid(image, vox_counts, True)
        #TODO construct centered grid
        output_scaled_mrcfile(ndimage.zoom(mrc.data, scale_factor), out_file, mrc)

def write_grid_to_mrcfil(grid, out_file, origin):
    """
    write out a scaled mrc file
        Args:
            grid (MRCImage): density at vertex grid
            out_file (pathlib.Path) the output file
            origin ([float]) the origin of the original grid
    """
    # TODO convert grid to np array & work out new origin
    # with mrcfile.new(out_file, overwrite=True) as mrc:
    #     mrc.set_data(data)
    #     mrc.header.origin = original_mrc.header.origin
    pass

def refine_volume_data(in_file, out_file, resolution):
    """
    simplify and smooth volume data in file
        Args:
            in_file (pathlib.Path) the input file
            out_file (pathlib.Path) the output file
            resolution (float) the new resolution
    """
    with mrcfile.open(in_file, mode='r+') as mrc:
        scale_factor=(mrc.voxel_size.tolist()[0])/resolution
        output_scaled_mrcfile(ndimage.zoom(mrc.data, scale_factor), out_file, mrc)

    return scale_factor

def output_scaled_mrcfile(data, out_file, original_mrc):
    """
    write out a scaled mrc file
        Args:
            data (numpy.ndarray): the volume data
            out_file (pathlib.Path) the output file
            original_mrc (mrcfile.File) the original input file
    """
    with mrcfile.new(out_file, overwrite=True) as mrc:
        mrc.set_data(data)
        mrc.header.origin = original_mrc.header.origin
        old_vox = original_mrc.voxel_size.tolist()
        new_list =[]
        nx_2 = mrc.header.nx
        ny_2 = mrc.header.ny
        nz_2 = mrc.header.nz
        scale_x = nx_2/original_mrc.header.nx
        scale_y = ny_2/original_mrc.header.ny
        scale_z = nz_2/original_mrc.header.nz
        new_list.append((1/scale_x)*old_vox[0]*(nx_2/(nx_2-1)))
        new_list.append((1/scale_y)*old_vox[1]*(ny_2/(ny_2-1)))
        new_list.append((1/scale_z)*old_vox[2]*(nz_2/(nz_2-1)))
        mrc.voxel_size = tuple(new_list)

        # not sure if doing different scale for each dimension necessary.
        # maybe fine to just use scale_x for all?
