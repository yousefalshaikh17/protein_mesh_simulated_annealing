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

      Vertex Indices:
         7+----------+6
         /|         /|
        / |        / |
      4+----------+5 |
       |  |       |  |         Axes:
       | 3+-------|--+2        z  y
       | /        | /          | /
       |/         |/           |/
      0+----------+1           +----x
"""

# set up linting
# pylint: disable = import-error

import sys
from time import process_time
from itertools import count
import numpy as np
import mrcfile

import ffeamesh.coord_utility as cu
from ffeamesh import utility
import ffeamesh.voxels2tets_utility as v2t
import ffeamesh.mrclininterpolate as mi

from ffeamesh.grid import Grid

def make_progress_test(end_x, end_y, end_z, steps=10, start_x=0, start_y=0, start_z=0):
    """
    make a progress test function
    Args:
        end_x (int): number of x columns
        end_y (int): number of y columns
        end_z (int): number of z columns
        steps (int): the number of iterations between reporting
        start_x (int): the count (zero start) of the first x column
        start_y (int): the count (zero start) of the first y column
        start_z (int): the count (zero start) of the first z column
    Returns:
        function(int, int, int) =>int or None
        (int): the total number of iterations
    """
    total = (end_x - start_x) * (end_y - start_y) * (end_z - start_z)
    stage = int(total/steps)

    def progress_test(current_iteration):
        """
        test if the current iteration should be reported for prograss
        Args:
            current_iteration (int): the number of the current iteration
        Returns:
            if the current iteration should be reported, the iteration , else
        """
        if current_iteration%stage == 0:
            print(f"Completed {current_iteration} out of {total} iterations")

    return progress_test

def all_voxels_to_tets(image, counts, progress, use_six_tets=True):
    """
    Converts image into voxels of 5 tetrohedrons.
    Args:
        image (MRCImage): the input file
        counts ([int, int, int]): voxel counts on x, y and z axis
        progress (bool): if true print out progress
        five_tets (bool): it True use 5 tet decomp, else 6
    Returns:
        (Grid)
    """
    # get the start and end value of the image cube axis
    start = [image.x_origin,
             image.y_origin,
             image.z_origin]
    end = [image.cell_size[0]+start[0],
           image.cell_size[1]+start[1],
           image.cell_size[2]+start[2]]
    image_counts =[image.get_nx(), image.get_ny(), image.get_nz()]

    grid = Grid(counts, start, end, image_counts)
    if progress:
        print("Grid object constructed.", file=sys.stdout)
    grid.build_grid(image)
    if progress:
        print("Vertices constructed", file=sys.stdout)
    if use_six_tets:
        grid.build_six_tets()
    else:
        grid.build_five_tets()
    if progress:
        print("Tets constructed", file=sys.stdout)

    return grid

def convert_mrc_to_tets(input_file,
                        output_file,
                        threshold,
                        ffea_out,
                        vtk_out,
                        verbose,
                        progress,
                        vox_counts,
                        low_vertices,
                        use_six_tets = False):
    """
    Converts the contents of a mrc file to a tetrohedron array (5 pre voxel).
    Args:
        input_file (pathlib.Path): name of input file
        output_file (pathlib.Path): name stem for output files
        threshold (float): the threshold below which results are ignored (isosurface value)
        ffea_out (bool): if true produce ffea input files (tetgen format)
        vtk_out (bool): if true produce vtk file
        verbose (bool): if true give details of results
        progress (bool): if true print out progress
        vox_counts ([int]): voxels counts on x, y, z, if None use image voxel counts
        low_vertices (int): number of vertices below isovalue that triggers culling
        use_six_tets (bool): if True convert grid voxels to six tets, else five
    Returns:
        None
    """
    prune_level = v2t.PruneLevel(low_vertices)

    # start time
    time_start = process_time()

    with mrcfile.mmap(input_file, mode='r+') as mrc:
        image = mi.MRCImage(mrc)

        if vox_counts is None:
            vox_counts = [image.get_nx(), image.get_ny(), image.get_nz()]

        grid = all_voxels_to_tets(image, vox_counts, progress, use_six_tets)
        if progress:
            print(f"Grid formed {len(grid.get_vertices())} vertices")

        grid.crop_mesh_to_isovalue(threshold, prune_level, progress)
        if progress:
            print(f"Grid formed {len(grid.get_connectivities())} tets")

        grid.remove_surplas_vertices()
        if progress:
            print(f"redundant vertices removed {len(grid.get_vertices())} remain")

        if grid.get_total_num_voxels() <= 0:
            print(f"Error: threshold value of {threshold} yielded no results", file=sys.stderr)
            sys.exit()

        # end time
        time_end = process_time()

        if progress:
            print(f"Conversion in {time_end - time_start} seconds, writing files")

        print(f"num voxels: {grid.get_total_num_voxels()}")
        print(f"num vertices {len(grid.get_vertices())}")
        print(f"num tets {len(grid.get_connectivities())}")

        if len(grid.get_connectivities()) == 0:
            print("No tetrahedrons were made, so no files written", file=sys.stdout)
            return

        if verbose:
            utility.verbose_output(mrc,
                                   grid.get_vertices(),
                                   grid.get_connectivities(),
                                   grid.get_total_num_voxels())

        v2t.write_tets_to_files(grid.get_vertices(),
                                grid.get_connectivities(),
                                output_file,
                                ffea_out,
                                vtk_out)
