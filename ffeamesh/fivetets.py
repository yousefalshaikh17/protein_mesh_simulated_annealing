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

import sys
import datetime
import getpass
from itertools import count
import numpy as np
import mrcfile
import vtk.util.numpy_support
import ffeamesh.coord_utility as cu
from ffeamesh.ffea_write import write_ffea_output
from ffeamesh import vtk_write
import ffeamesh.vtk_utility as vtk_u
from ffeamesh import utility

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

def voxel_into_5_tets(voxel_x, voxel_y, voxel_z, frac_to_cart, coord_store, connectivities_final):
    """
    convert a single voxel into 5 tets
    Args:

    """
    coords = []
    create_cube_coords(voxel_x, voxel_y, voxel_z, frac_to_cart, coords)

    indices = None
    if is_odd(voxel_x, voxel_y, voxel_z):
        indices = odd_cube_tet_indecies()
    else:
        indices = even_cube_tet_indices()

        # test the tets and append those that pass
    for tet in indices:
        tet_indices = []
        for index in tet:
            tet_indices.append(coord_store.add(coords[index]))
        connectivities_final.append(tet_indices)

def voxels_to_5_tets_plain(mrc, threshold, progress):
    """
    converts voxels to tetrohedrons and returns those with average values above a threshold
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
    frac_to_cart = make_fractional_to_cartesian_conversion_function(mrc)
    prog_test = make_progress_test(mrc.header.nx, mrc.header.ny, mrc.header.nz)

    # Create an array of array of 8 point (co-ordinates) for each hexahedron (voxel)
    for voxel_z in range(0, mrc.header.nz):
        for voxel_y in range(0, mrc.header.ny):
            for voxel_x in range(0, mrc.header.nx):
                if progress:
                    prog_test(next(curent_iteration))

                # Threshold the voxels out of the mrc map data
                if mrc.data[voxel_z, voxel_y, voxel_x] > threshold:
                    voxel_into_5_tets(voxel_x, voxel_y, voxel_z, frac_to_cart, coord_store, connectivities_final)
                    next(voxel_count)

    points = [[coord.cart.x, coord.cart.y, coord.cart.z] for coord in coord_store.to_list()]

    return next(voxel_count), points, connectivities_final

def convert_mrc_to_5tets(input_file, output_file, threshold, ffea_out, vtk_out, verbose, progress):
    """
    Converts the contents of an mrc file to a tetrohedron array.
    It is called by the fivetets.py script and controls the whole conversion.
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
            print(f"number of voxels over threshold {nvoxel}")
            utility.print_voxel_stats(points, connectivities_final)

        write_tets_to_files(points, connectivities_final, output_file, ffea_out, vtk_out)

def convert_mrc_to_5tets_interp(input_file, output_file, threshold, ffea_out, vtk_out, verbose, progress):
    """
    Converts the contents of a mrc file to a tetrohedron array.
    Is called from the five_tets.py script and controls the whole conversion.
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
        nvoxel, points, tet_connectivities, = voxels_to_5_tets_interp(mrc, threshold, progress)

        if nvoxel <= 0:
            print(f"Error: threshold value of {threshold} yielded no results", file=sys.stderr)
            sys.exit()

        if verbose:
            print(f"number of voxels over threshold {nvoxel}")
            utility.print_voxel_stats(points, tet_connectivities)

        write_tets_to_files(points, tet_connectivities, output_file, ffea_out, vtk_out)

def voxels_to_5_tets_interp(mrc, threshold, progress):
    """
    converts voxels to tetrohedrons and returns those with average values above a threshold
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
    frac_to_cart = make_fractional_to_cartesian_conversion_function(mrc)
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
                cube_vertex_values = make_vertex_values(voxel_x, voxel_y, voxel_z, mrc)

                # test is at least one is over the the threshold
                if sum([np.count_nonzero(x>threshold) for x in cube_vertex_values]) > 0:
                    # count the number of voxels
                    next(voxel_count)

                    # make the matching coordinates
                    coords = []
                    create_cube_coords(voxel_x, voxel_y, voxel_z, frac_to_cart, coords)

                    # connectivity of 5 tets in single voxel
                    indices = None
                    if is_odd(voxel_x, voxel_y, voxel_z):
                        indices = odd_cube_tet_indecies()
                    else:
                        indices = even_cube_tet_indices()

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
                            connectivities_final.append(tet_indices)

    tmp = [[coord.cart.x, coord.cart.y, coord.cart.z] for coord in coord_store.to_list()]
    return next(voxel_count), tmp, connectivities_final

def even_cube_tet_indices():
    """
    return a list of lists for the constuction of 5 tets from the 8 vertices of an even cube
    Args:
        None
    Returns
        lits(list) : five lists of four vertex indices representing the tets
    """
    return [[0, 4, 5, 7],
            [0, 1, 2, 5],
            [2, 5, 6, 7],
            [0, 2, 3, 7],
            [0, 2, 5, 7]]

def even_cube_tets(cube):
    """
    convert a list of the eight vertices of a an even cube
    into a list of lists representing tets
    Args:
        cube (list): the eight vertices of a cube
    Returns
        lits(list) : five lists of four vertices representing the tets
    """
    tet1 = np.array([cube[0], cube[4], cube[5], cube[7]])
    tet2 = np.array([cube[0], cube[1], cube[2], cube[5]])
    tet3 = np.array([cube[2], cube[5], cube[6], cube[7]])
    tet4 = np.array([cube[0], cube[2], cube[3], cube[7]])
    tet5 = np.array([cube[0], cube[2], cube[5], cube[7]])

    tet_list = [tet1, tet2, tet3, tet4, tet5]

    return tet_list

def odd_cube_tet_indecies():
    """
    return a list of lists for the constuction of 5 tets from the 8 vertices of an odd cube
    Args:
        None
    Returns
        lits(list) : five lists of four vertex indices representing the tets
    """
    return [[0, 1, 3, 4],
            [1, 4, 5, 6],
            [1, 2, 3, 6],
            [3, 4, 6, 7],
            [1, 3, 4, 6]]

def odd_cube_tets(cube):
    """
    convert a list of the eight vertices of a an odd cube
    into a list of lists representing tets
    Args:
        cube (list): the eight vertices of a cube
    Returns
        lits(list) : five lists of four vertices representing the tets
    """
    tet1 = np.array([cube[0], cube[1], cube[3], cube[4]])
    tet2 = np.array([cube[1], cube[4], cube[5], cube[6]])
    tet3 = np.array([cube[1], cube[2], cube[3], cube[6]])
    tet4 = np.array([cube[3], cube[4], cube[6], cube[7]])
    tet5 = np.array([cube[1], cube[3], cube[4], cube[6]])

    tet_list = [tet1, tet2, tet3, tet4, tet5]

    return tet_list

def is_odd(x_index, y_index, z_index):
    '''
    Logic to decide if a voxel is an odd or and even.
    Args:
        x_index (int) Index for x value in values array.
        y_index (int) Index for y value in values array.
        z_index (int) Index for z value in values array.
    Returns:
        (bool)  True if voxel in an odd position, else odd
    '''
    flag = None
    y_parity = (y_index % 2 == 0)
    x_parity = (x_index % 2 == 0)

    # Logic for alternating tet division 0 (even) or 1 (odd) - to identify the odd and even voxels
    if z_index % 2 == 0:
        if y_parity:
            flag = not x_parity
        else:
            flag = x_parity
    else:
        if y_parity:
            flag = x_parity
        else:
            flag = not x_parity

    return flag

def make_fractional_to_cartesian_conversion_function(mrc):
    """
    make_fractional_to_cartesian_conversion_function(mrc)
    make a functor that will convert a point defined by fractional
    coordinates in the unit cell to cartesian coordinates
    TODO imp for angles not 90 degrees
    Args:
        mrc (mrcfile): the mrc file
    Returns
        function: fractional (i, j, k) (flaot) => cartesian (x, y, z) (float)
    """
    # get cell sizes
    x_res = mrc.header.cella.x/np.float32(mrc.header.mx)
    y_res = mrc.header.cella.y/np.float32(mrc.header.my)
    z_res = mrc.header.cella.z/np.float32(mrc.header.mz)

    # get coordinate origine
    x_trans = np.float32(mrc.header.origin['x'])
    y_trans = np.float32(mrc.header.origin['y'])
    z_trans = np.float32(mrc.header.origin['z'])

    def fractional_to_cartesian_coordinates(x_frac, y_frac, z_frac):
        """
        convert integer indices to 3D coordinates.
        Args:
            x_frac (float): the fractional x coordinate
            y_frac (float): the fractional y coordinate
            z_frac (float): the fractional z coordinate
        Returns:
            (CoordTransform): container for fractional and cartesian coordinates
        """
        x_cart = (x_frac * x_res) + x_trans
        y_cart = (y_frac * y_res) + y_trans
        z_cart = (z_frac * z_res) + z_trans

        return cu.CoordTransform(cu.Coordinate(x_frac, y_frac, z_frac),
                              cu.Coordinate(x_cart, y_cart, z_cart))

    return fractional_to_cartesian_coordinates

def create_cube_coords(x_index, y_index, z_index, frac_to_cart, coords):
    '''
    Caluculates the next 8 coords for the next volxel that has been thresholded
    previously (logic in loop that calls this one).
    Args:
        x_index (int):          Indecies to the x coord in the coords array.
        y_index (int)           Indecies to the y coord in the coords array.
        z_index (int)           Indecies to the z coord in the coords array.
        frac_to_cart (float*3=>CoordTransform): fractional to cartesian conversin function
        coords (CoordTransform list): The coordinates of the thresholded voxels passed by
                                      reference to this function.

    Returns:
        None
    '''
    # Calculate the cordinates of each vertex of the voxel
    coords.append(frac_to_cart((x_index-0.5), (y_index-0.5), (z_index-0.5)))
    coords.append(frac_to_cart((x_index+0.5), (y_index-0.5), (z_index-0.5)))
    coords.append(frac_to_cart((x_index+0.5), (y_index+0.5), (z_index-0.5)))
    coords.append(frac_to_cart((x_index-0.5), (y_index+0.5), (z_index-0.5)))
    coords.append(frac_to_cart((x_index-0.5), (y_index-0.5), (z_index+0.5)))
    coords.append(frac_to_cart((x_index+0.5), (y_index-0.5), (z_index+0.5)))
    coords.append(frac_to_cart((x_index+0.5), (y_index+0.5), (z_index+0.5)))
    coords.append(frac_to_cart((x_index-0.5), (y_index+0.5), (z_index+0.5)))

def make_vertex_values(x_index, y_index, z_index, mrc):
    """
    make an arry of eight values by averaging at the vertices
    Args:
        x_index (int)       Indecies to the x coord in the coords array.
        y_index (int)       Indecies to the y coord in the coords array.
        z_index (int)       Indecies to the z coord in the coords array.
        mrc (mrcfile.mmap)  image source file
    Returns:
        (float list): eight values
    """
    ctr_value = mrc.data[z_index, y_index, x_index]
    values = []

    # compute the eight averages
    values.append((ctr_value + mrc.data[z_index-1, y_index-1, x_index-1])/2.0)
    values.append((ctr_value + mrc.data[z_index-1, y_index-1, x_index+1])/2.0)
    values.append((ctr_value + mrc.data[z_index-1, y_index+1, x_index+1])/2.0)
    values.append((ctr_value + mrc.data[z_index-1, y_index+1, x_index-1])/2.0)
    values.append((ctr_value + mrc.data[z_index+1, y_index-1, x_index-1])/2.0)
    values.append((ctr_value + mrc.data[z_index+1, y_index-1, x_index+1])/2.0)
    values.append((ctr_value + mrc.data[z_index+1, y_index+1, x_index+1])/2.0)
    values.append((ctr_value + mrc.data[z_index+1, y_index+1, x_index+1])/2.0)

    return values

def write_tets_to_files(points_list, tets_connectivity, output_file, ffea_out, vtk_out):
    """
    Sets it up so data goes to the ffea writer and the vtk writer correctly.
    Args:
        points_list (float*3 list): 3d coordinates of the points forming the tets
        tets_connectivity (int*4 list): for each tet the indices of its vertices in the points list
        output_file (pathlib.Path): name stem for ouput files
        ffea_out (bool): if true write ffea input files
        vtk_out (bool): if true write a vtk file
    Returns:
        None
    """
    # convert coords to np.array
    points_np = np.array(points_list)

    cells_con = vtk_u.make_vtk_tet_connectivity(tets_connectivity)

    # make the grid (vtk scene)
    vtk_pts = vtk.vtkPoints()
    vtk_pts.SetData(vtk.util.numpy_support.numpy_to_vtk(points_np, deep=True))
    grid = vtk.vtkUnstructuredGrid() #create unstructured grid
    grid.SetPoints(vtk_pts) #assign points to grid
    grid.SetCells(vtk.VTK_TETRA, cells_con) #assign tet cells to grid

    #write vtk file
    if vtk_out:
        vtk_write.vtk_output(grid, output_file)

    # write tetgen file for ffea input
    if ffea_out:
        ffea_output(grid, points_np, tets_connectivity, output_file)

def ffea_output(grid, points, tets_connectivity, output_file):
    """
    Constructs the faces and calls the ffea writer to output file.
    Args:
        grid (vtk.vtkUnstructuredGrid): vtk scene
        points (float np.ndarray): duplicate free list of vertices
        tets_connectivity (int*4 list): for each tet the indices of its vertices in the points list
        output_file (pathlib.Path): name stem of output files
    Returns:
        None
    """
    surf_points, cell_data, cell_count = vtk_u.get_vtk_surface(grid)

    # make a connectivity into the tets points
    original_ids = np.zeros((len(surf_points),), dtype='int16')
    for pos, point in enumerate(surf_points):
        # index of surface point in the surface faces's points array
        original_ids[pos] = np.where((points==point).all(axis=1))[0]

    # This holds true if all polys are of the same kind, e.g. triangles.
    assert cell_data.GetNumberOfValues()%cell_count==0

    # reshape the cells array to match tet gen output standard
    col_count = cell_data.GetNumberOfValues()//cell_count
    numpy_cells = np.array(cell_data)
    faces = numpy_cells.reshape((-1, col_count))

    #write to tetgen .ele, .node, .face
    date = datetime.datetime.now().strftime("%x")
    write_ffea_output(output_file,
                      tets_connectivity,
                      points,
                      faces,
                      original_ids,
                      f'# created by {getpass.getuser()} on {date}')
