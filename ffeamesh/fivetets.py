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
from itertools import count
import numpy as np
import mrcfile
import ffeamesh.coord_utility as cu
from ffeamesh import utility
import ffeamesh.voxels2tets_utility as v2t
import ffeamesh.mrclininterpolate as mi

def convert_mrc_to_5tets(input_file, output_file, threshold, ffea_out, vtk_out, verbose, progress):
    """
    Converts the contents of an mrc file to a tetrohedron array (5 tets pre voxel).
    Args:
        input_file (pathlib.Path): name of input file
        output_file (pathlib.Path): name stem for output files
        threshold (float): the threshold below which results are ignored (isosurface value)
        ffea_out (bool): if true produce ffea input files (tetgen format)
        vtk_out (bool): if true produce vtk file
        verbose (bool): if true give details of results
        progress (bool): if true print out progress
    Returns:
        None
    """
    # Reads mrc file into a map which is the fastest way to work with the data and has
    # a header and data section in it.
    with mrcfile.mmap(input_file, mode='r+') as mrc:
        nvoxel, points, connectivities_final = voxels_to_5_tets_plain(mrc, threshold, progress)

        if nvoxel < 1:
            print(f"Error: threshold value of {threshold} yielded no voxels", file=sys.stderr)
            sys.exit()

        if verbose:
            utility.verbose_output(mrc, points, connectivities_final, nvoxel)

        v2t.write_tets_to_files(points, connectivities_final, output_file, ffea_out, vtk_out)

def convert_mrc_to_5tets_interp(input_file,
                                output_file,
                                threshold,
                                ffea_out,
                                vtk_out,
                                verbose,
                                progress):
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
    Returns:
        None
    """
    with mrcfile.mmap(input_file, mode='r+') as mrc:
        nvoxel, points, tet_connectivities = voxels_to_5_tets_interp(mrc, threshold, progress)

        if nvoxel <= 0:
            print(f"Error: threshold value of {threshold} yielded no results", file=sys.stderr)
            sys.exit()

        if verbose:
            utility.verbose_output(mrc, points, tet_connectivities, nvoxel)

        v2t.write_tets_to_files(points, tet_connectivities, output_file, ffea_out, vtk_out)

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

def plain_voxel_to_5_tets(voxel, frac_to_cart, coord_store, connectivities):
    """
    convert a single voxel into 5 tets
    Args:
        voxel (Coordinate): array indices of voxel
        frac_to_cart (function fractional => cartesian): coordinages transform
        coord_store (CoordTransform, []): list for coordinates
        connectivities ([[int, int, int int]]): tet connectivities
    Returns:
        None
    """
    coords = v2t.create_cube_coords(voxel, frac_to_cart)

    indices = None
    if v2t.is_odd(voxel.x, voxel.y, voxel.z):
        indices = v2t.odd_cube_tet_indices()
    else:
        indices = v2t.even_cube_tet_indices()

        # test the tets and append those that pass
    for tet in indices:
        tet_indices = []
        for index in tet:
            tet_indices.append(coord_store.add(coords[index]))
        connectivities.append(tet_indices)

def voxels_to_5_tets_plain(mrc, threshold, progress):
    """
    Thresholds voxels and converts them all into 5 tetrohedrons.
    Args:
        mrc (mrcfile.mmap): the input file
        threshold (float): the acceptance limit
        progress (bool): if true print out progress
    Returns:
        (int): the number of voxels transformed
        ([float, float, float] list): the coordinates of the vertices
        ([int, int, int int] list): the vertices of the tets as indices in the coordinates list
    """
    coord_store = cu.UniqueTransformStore()
    connectivities_final = []
    voxel_count = count(0)
    curent_iteration = count(1)
    frac_to_cart = v2t.make_fractional_to_cartesian_conversion_function(mrc)
    prog_test = make_progress_test(mrc.header.nx, mrc.header.ny, mrc.header.nz)

    # Create an array of array of 8 point (co-ordinates) for each hexahedron (voxel)
    for voxel_z in range(0, mrc.header.nz):
        for voxel_y in range(0, mrc.header.ny):
            for voxel_x in range(0, mrc.header.nx):
                if progress:
                    prog_test(next(curent_iteration))

                # Threshold the voxels out of the mrc map data
                if mrc.data[voxel_z, voxel_y, voxel_x] > threshold:
                    plain_voxel_to_5_tets(cu.Coordinate(voxel_x, voxel_y, voxel_z),
                                          frac_to_cart,
                                          coord_store,
                                          connectivities_final)
                    next(voxel_count)

    points = [[coord.cart.x, coord.cart.y, coord.cart.z] for coord in coord_store.to_list()]

    return next(voxel_count), points, connectivities_final

def voxels_to_5_tets_interp(mrc, threshold, progress):
    """
    Converts voxels into 5 tetrohedrons and returns only the tetrahedrons above a threshold.
    Args:
        mrc (mrcfile.mmap): the input file
        threshold (float): the acceptance limit
        progress (bool): if true print out progress
    Returns:
        (int): the number of voxels transformed
        ([float, float, float] list): the coordinates of the vertices
        ([int, int, int int] list): the vertices of the tets as indices in the coordinates list
    """
    coord_store = cu.UniqueTransformStore()
    connectivities_final = []
    curent_iteration = count(1)
    voxel_count = count(0)
    frac_to_cart = v2t.make_fractional_to_cartesian_conversion_function(mrc)
    prog_test = make_progress_test(mrc.header.nx-1,
                                   mrc.header.ny-1,
                                   mrc.header.nz-1,
                                   start_x=1,
                                   start_y=1,
                                   start_z=1)

    for voxel_z in range(1, mrc.header.nz-1):
        for voxel_y in range(1, mrc.header.ny-1):
            for voxel_x in range(1, mrc.header.nx-1):
                if progress:
                    prog_test(next(curent_iteration))

                # find the interpolated values at the vertices
                cube_vertex_values = v2t.make_vertex_values(voxel_x, voxel_y, voxel_z, mrc)

                # test is at least one is over the the threshold
                if sum(np.count_nonzero(x>threshold) for x in cube_vertex_values) > 0:
                    # count the number of voxels
                    next(voxel_count)

                    interp_voxel_to_5_tets(cu.Coordinate(voxel_x, voxel_y, voxel_z),
                                           frac_to_cart,
                                           threshold,
                                           cube_vertex_values,
                                           coord_store,
                                           connectivities_final)

    tmp = [[coord.cart.x, coord.cart.y, coord.cart.z] for coord in coord_store.to_list()]
    return next(voxel_count), tmp, connectivities_final

def interp_voxel_to_5_tets(voxel,
                           frac_to_cart,
                           threshold,
                           cube_vertex_values,
                           coord_store,
                           connectivities):
    """
    Convert a single voxel to 5 tets thresholding the tets by interpolating
    vertex values. Called by voxels_to_6_tets_interp
    Args:
        voxel (Coordinate): the array indices of the voxel
        frac_to_cart (fractional => cartesian): coordinates conversion function
        threshold (float): threshold for inclusion
        cube_vertex_values ([float]): the image values interpolated at voxel vertices
        coord_store (UniqueTransformStore): store for the vertices
        connectivities ([[int, int, int, int]]): the tet's connectivites
    Returns:
        None
    """
    # make the matching coordinates
    coords = v2t.create_cube_coords(voxel, frac_to_cart)

    # connectivity of 5 tets in single voxel
    indices = None
    if v2t.is_odd(voxel.x, voxel.y, voxel.z):
        indices = v2t.odd_cube_tet_indices()
    else:
        indices = v2t.even_cube_tet_indices()

    # test the tets and append those that pass
    for tet in indices:
        average = 0.0
        for index in tet:
            average += cube_vertex_values[index]
        average /= 4
        if average > threshold:
            tet_indices = []
            for index in tet:
                tet_indices.append(coord_store.add(coords[index]))
            connectivities.append(tet_indices)

def all_voxels_to_5_tets(mrc, counts, progress):
    """
    Converts image into voxels of 5 tetrohedrons.
    Args:
        mrc (mrcfile.mmap): the input file
        counts ([int, int, int]): voxel counts on x, y and z axis
        progress (bool): if true print out progress
    Returns:
        (int): the number of voxels transformed
        ([float, float, float] list): the coordinates of the vertices
        ([int, int, int int] list): the vertices of the tets as indices in the coordinates list
    """
    x_lin = np.linspace(0.0, mrc.header.cella.x, counts[0])
    y_lin = np.linspace(0.0, mrc.header.cella.y, counts[1])
    z_lin = np.linspace(0.0, mrc.header.cella.z, counts[2])
    print(x_lin, y_lin, z_lin)
    return None, None, None

def convert_mrc_to_5tets_interp2(input_file,
                                output_file,
                                threshold,
                                ffea_out,
                                vtk_out,
                                verbose,
                                progress,
                                vox_counts):
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
        vox_counts ([int]): numbers of voxels on each dimension x, y, z, if None use image
    Returns:
        None
    """
    with mrcfile.mmap(input_file, mode='r+') as mrc:
        nvoxel=None
        points=None
        tet_connectivities=None
        if vox_counts is None:
            nvoxel, points, tet_connectivities = voxels_to_5_tets_interp(mrc, threshold/2.0, progress)
        else:
            nvoxel, points, tet_connectivities = all_voxels_to_5_tets(mrc, vox_counts, progress)
            quit(f"by now: {nvoxel} {points} {tet_connectivities}")

        if nvoxel <= 0:
            print(f"Error: threshold value of {threshold} yielded no results", file=sys.stderr)
            sys.exit()

        image = mi.MRCImage(mrc)
        connectivity = v2t.crop_mesh_to_isovalue(points, tet_connectivities, image, threshold)

        if verbose:
            utility.verbose_output(mrc, points, connectivity, nvoxel)

        v2t.write_tets_to_files(points, connectivity, output_file, ffea_out, vtk_out)
