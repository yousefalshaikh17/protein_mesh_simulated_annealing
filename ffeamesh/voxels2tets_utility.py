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
import datetime
import getpass
import numpy as np
import enum
import vtk
from vtkmodules.util import numpy_support
from ffeamesh import vtk_write
import ffeamesh.coord_utility as cu
from ffeamesh.ffea_write import write_ffea_output
import ffeamesh.vtk_utility as vtk_u

class PruneLevel(enum.IntEnum):
    """
    number of verticies allowed below isovlue
    """
    FOUR  = 4
    THREE = 3
    TWO   = 2
    ONE   = 1

def make_fractional_to_cartesian_conversion_function(mrc):
    """
    make_fractional_to_cartesian_conversion_function(mrc)
    make a functor that will convert a point defined by fractional
    coordinates in the unit cell to cartesian coordinates
    TODO imp for angles not 90 degrees
    Args:
        mrc (mrcfile): the mrc file
    Returns
        function: fractional (i, j, k) (float) => cartesian (x, y, z) (float)
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

def cube_6_tet_indices():
    """
    return a list of lists for the constuction of 6 tets from the 8 vertices a cube
    Args:
        None
    Returns
        lits(list) : five lists of four vertex indices representing the tets
    """
    return [[0, 6, 5, 1],
            [0, 6, 1, 2],
            [0, 6, 2, 3],
            [0, 6, 3, 7],
            [0, 6, 7, 4],
            [0, 6, 4, 5]]

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
    indices_array = even_cube_tet_indices()
    tet_list = []
    for indices in indices_array:
        tet_list.append(np.array([cube[indices[0]],
                                  cube[indices[1]],
                                  cube[indices[2]],
                                  cube[indices[3]]]))

    return tet_list

def odd_cube_tet_indices():
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
    indices_array = odd_cube_tet_indices()
    tet_list = []
    for indices in indices_array:
        tet_list.append(np.array([cube[indices[0]],
                                  cube[indices[1]],
                                  cube[indices[2]],
                                  cube[indices[3]]]))

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
    y_parity =  y_index % 2 == 0
    x_parity =  x_index % 2 == 0

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

def create_cube_coords(voxel, frac_to_cart):
    '''
    Caluculates the next 8 coords for the next voxel that has been
    previously thresholded (logic in loop that calls this one).
    Args:
        voxel (Coordinate): array indices of voxel
        frac_to_cart (float*3=>CoordTransform): fractional to cartesian conversin function
    Returns:
        None
    '''
    coords = []
    # Calculate the cordinates of each vertex of the voxel
    coords.append(frac_to_cart((voxel.x-0.5), (voxel.y-0.5), (voxel.z-0.5)))
    coords.append(frac_to_cart((voxel.x+0.5), (voxel.y-0.5), (voxel.z-0.5)))
    coords.append(frac_to_cart((voxel.x+0.5), (voxel.y+0.5), (voxel.z-0.5)))
    coords.append(frac_to_cart((voxel.x-0.5), (voxel.y+0.5), (voxel.z-0.5)))
    coords.append(frac_to_cart((voxel.x-0.5), (voxel.y-0.5), (voxel.z+0.5)))
    coords.append(frac_to_cart((voxel.x+0.5), (voxel.y-0.5), (voxel.z+0.5)))
    coords.append(frac_to_cart((voxel.x+0.5), (voxel.y+0.5), (voxel.z+0.5)))
    coords.append(frac_to_cart((voxel.x-0.5), (voxel.y+0.5), (voxel.z+0.5)))

    return coords

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

    # connectivity in vtk
    cells_con = vtk_u.make_vtk_tet_connectivity(tets_connectivity)

    # make the grid (vtk scene)
    vtk_pts = vtk.vtkPoints()

    # work on deep copy
    vtk_pts.SetData(numpy_support.numpy_to_vtk(points_np, deep=True))

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

###########################################################################
###########################################################################

def find_indices_of(point, points_array):
    """
    find all indices of point in array
    Args:
        point [float, float, float]: target point
        points_array [[float, float, float]...]: array of points
    Retuns
        [int]: list of indices of occurance of point in points_array
    """
    indices = []
    for index, pt in enumerate(points_array):
        if point[0] == pt[0] and point[1] == pt[1] and point[2] == pt[2]:
            indices.append(index)

    if len(indices) == 0:
        raise ValueError(f"point not in list {point}")

    return indices

def crop_mesh_to_isovalue(points_list, tets_connectivity, image, isovalue):
    """
    Remove tets outside or largly outside isovalue
    Args:
        points_list (float*3 list): 3d coordinates of the points forming the tets
        tets_connectivity (int*4 list): for each tet the indices of its vertices in the points list
        image (MRCImage): image for interpolation
        isovalue (float): limit value
    Returns:
        TODO
    """
    new_tets = find_tets_outside_isosurface(points_list, tets_connectivity, image, isovalue)

    print(f"Number of new tets is {len(new_tets)}\n######################\n")

    return new_tets


def find_tets_outside_isosurface(points_list, tets_connectivity, image, isovalue):
    """
    Remove surface tets outside or with one face outside isovalue
    Args:
        points_list (float*3 list): 3d coordinates of the points forming the tets
        tets_connectivity (int*4 list): for each tet the indices of its vertices in the points list
        image (MRCImage): image for interpolation
        isovalue (float): limit value
    Returns:
        (int*4 list): for each tet passing test the indices of its vertices in the points list
    """
    # convert coords to np.array
    points_np = np.array(points_list)

    return prune_mesh(points_np, tets_connectivity, image, isovalue)

def make_vtk_grid(points_np, tets_connectivity):
    """
    make a vtk UnstructuredGrid
    Args:
        points_list (float*3 numpy.array): 3d coordinates of the points forming the tets
        tets_connectivity (int*4 list): for each tet the indices of its vertices in the points list
    Returns:
        UnstructuredGrid
    """
    # convert coords to vtk.vtkCellArray
    cells_con = vtk_u.make_vtk_tet_connectivity(tets_connectivity)

    # make contaner fo the points
    vtk_pts = vtk.vtkPoints()

    # copy the numpy points into the vtk container
    vtk_pts.SetData(numpy_support.numpy_to_vtk(points_np, deep=True))

    #create unstructured grid (conatiner for any combinations of any cell types)
    grid = vtk.vtkUnstructuredGrid()
    #assign points to grid
    grid.SetPoints(vtk_pts)
    #assign tet cells to grid
    grid.SetCells(vtk.VTK_TETRA, cells_con)

    return grid

def crop_mesh_to_isovalue(vertices,  densities, connectivity, isovalue, level=PruneLevel.TWO):
    """
    Remove tets outside or largly outside isovalue
    Args:
        vertices ([[float, float float]]): 3d coordinates of the points forming the tets
        densities (float): the image density at each vertex
        connectivity ([[int, int, int, int]]): the indices of each tet's vertices in the list
        isovalue (float): limit value
        level (PruneLevel): the number of vertices below the isovalue that causes deleation
    Returns:
        [[int, int, int, int]]: connectivity of tets passing test
    """
    print(f"prune to {isovalue}")
    print(f"Number of points {len(vertices)}")
    print(f"Number of tets {len(connectivity)}")
    print(f"Prune level {level}")

    # find tets to be removed
    ##########################################################
    tets_for_deletion = []
    out_count = 0
    # out_4 = 0
    # out_3 = 0
    # out_2 = 0
    # out_1 = 0

    for tet_index, tet in enumerate(connectivity):
        count_outside = 0
        for index in tet:
            if densities[index] < isovalue:
                count_outside += 1

        if count_outside > level.value:
            tets_for_deletion.append(tet_index)
            out_count += 1

        # if count_outside > 3:
        #     tets_for_deletion.append(tet_index)
        #     out_4 += 1
        # elif count_outside > 2:
        #     tets_for_deletion.append(tet_index)
        #     out_3 += 1
        # elif count_outside > 1:
        #     tets_for_deletion.append(tet_index)
        #     out_2 += 1
        # elif count_outside > 0:
        #     tets_for_deletion.append(tet_index)
        #     out_1 += 1

    # print(f"4 out {out_4}: 3 out {out_3}: 2 out {out_2}: 1 out {out_1}")
    print(f"Out count: {out_count}")

    # make new tet connectivity list
    ##########################################################
    new_tets = []
    for index, tet in enumerate(connectivity):
        if index not in tets_for_deletion:
            new_tets.append(tet)

    print("returning new tets")
    return new_tets
