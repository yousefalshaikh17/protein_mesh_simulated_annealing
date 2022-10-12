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
from collections import (namedtuple, OrderedDict)
import numpy as np
import mrcfile
import vtk.util.numpy_support
from ffeamesh.writers import write_ffea_output

## data structure for a 3D coordinate
Coordinate = namedtuple("Coordinate", "x, y, z")

class CoordTransform(object):
    """
    store for the fractional and cartesian coordinates of a 3D lattice point
    """
    def __init__(self, frac, cart):
        """
        initalize the object
        Args:
            frac (Coordinate): fractional coordinate
            cart (Coordinate): cartesian coordinate
        """
        ## storage for the fractional coordinate
        self.frac = frac

        ## storage for the cartesian coordinate
        self.cart = cart

    def __eq__(self, lhs):
        """
        equality is based on fractional coordinate only
        Args:
            lhs (CoordTransform): the object for comparison
        Returns:
            (bool) True if fractional coordinates are equal else False
        """
        if self.frac.x != lhs.frac.x:
            return False

        if self.frac.y != lhs.frac.y:
            return False

        if self.frac.z != self.frac.z:
            return False

        return True

    def __hash__(self):
        """
        hash value is based on the fractional coordinate only
        Returns:
            (int) the hash value of the fractional coordinate
        """
        return hash((self.frac.x, self.frac.y, self.frac.z))

    def __repr__(self):
        """
        string representation of the object
        """
        sf = f"({self.frac.x}, {self.frac.y})"
        sc = f"({self.cart.x}, {self.cart.y})"
        return f"<CoordTransform: {sf} => {sc}>"

class UniqueTransformStore(object):
    """
    store for CoordTransform which only holds unique objects,
    based on fractional coordinate
    """
    def __init__(self):
        """
        initalize the object
        """
        ## current size of the array
        self.current_size = 0

        ## dictionary holding the data, will work python < 3.7
        self.data = OrderedDict()

    def add(self, coord):
        """
        add a new CoordTransform
        Args:
            coord (CoordTransform): the object to be added
        Returns
            (int): the index of the object in the list of unique objects
        """
        if coord in self.data.keys():
            return self.data[coord]

        self.data[coord] = self.current_size
        tmp = self.current_size
        self.current_size += 1

        return tmp

    def to_list(self):
        """
        convert the keys to a list in entry order
        Returns:
            (CoordTransform list)
        """
        return [coord for coord in self.data.keys()]

    def __str__(self):
        """
        provide string representation
        """
        return f"<UniqueArray: {self.current_size} items>"

    def __hash__(self):
        """
        hash function is the hash of the data
        """
        return hash(self.data)

def convert_mrc_to_5tets_interp(input_file, output_file, threshold, ffea_out, vtk_out):
    """
    Converts the contents of an mrc file to a tetrohedron array.
    Is called from the five_tets.py script and controls the whole conversion.
    Args:
        input_file (pathlib.Path)
        output_file (pathlib.Path)
        threshold (float)
        ffea_out (bool)
        vtk_out (bool)
    Returns:
        None
    """
    with mrcfile.mmap(input_file, mode='r+') as mrc:
        nvoxel, coords_final, connectivities_final, = voxels_to_5_tets_thershold(mrc, threshold)

        if nvoxel <= 0:
            print(f"Error: threshold value of {threshold} yielded no voxels", file=sys.stderr)
            sys.exit()

        print(f"number of voxels over thershold {nvoxel}")

        cells_con = make_vtk_tet_connectivity(connectivities_final)

        # make the grid (vtk scene)
        vtk_pts = vtk.vtkPoints()
        vtk_pts.SetData(vtk.util.numpy_support.numpy_to_vtk(coords_final, deep=True))
        grid = vtk.vtkUnstructuredGrid() #create unstructured grid
        grid.SetPoints(vtk_pts) #assign points to grid
        grid.SetCells(vtk.VTK_TETRA, cells_con) #assign tet cells to grid

        #write vtk file
        vtk_output(grid, output_file)

        # make the connectivity data for voxels
        # points, cells = make_voxel_connectivity(nvoxel, coords)

        # make the connectivity for tets
        #tet_array = make_tet_connectivity(nvoxel, cells, alternate)

        #write_tets_to_files(nvoxel, points, cells, tet_array, output_file, ffea_out, vtk_out)

def voxels_to_5_tets_thershold(mrc, threshold):
    """
    """
    voxel_count = count(0)
    frac_to_cart = make_fractional_to_cartesian_conversion_function(mrc)

    coords_final = []
    connectivities_final = []

    for voxel_z in range(1, mrc.header.nz-1):
        for voxel_y in range(1, mrc.header.ny-1):
            for voxel_x in range(1, mrc.header.nx-1):
                # find the interpolated values at the vertices
                cube_vertex_values = make_vertex_values(voxel_x, voxel_y, voxel_z, mrc)

                # test is at least one is over the the threshold
                if sum([np.count_nonzero(x>threshold) for x in cube_vertex_values]) > 0:
                    # count the number of voxels
                    next(voxel_count)

                    # make the matching coordinates
                    coords = []
                    create_cube_coords(voxel_x, voxel_y, voxel_z, frac_to_cart, coords)

                    coord_end_index = len(coords_final)
                    coords_final += (coords)

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
                            connectivities_final.append(
                                [connect + coord_end_index for connect in tet])

    return next(voxel_count), coords_final, connectivities_final

def even_cube_tet_indices():
    """
    return a list of lists for the constuctio of 5 tets from the 8 vertices of an even cube
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
    return a list of lists for the constuctio of 5 tets from the 8 vertices of an odd cube
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
    # Logic for alternating tet division 0 (even) or 1 (odd) - to identify the odd and even voxels
    if z_index % 2 == 0:
        if y_index % 2 == 0:
            if x_index % 2 == 0:
                flag = False
            else:
                flag = True
        else:
            if x_index % 2 == 0:
                flag = True
            else:
                flag = False
    else:
        if y_index % 2 == 0:
            if x_index % 2 == 0:
                flag = True
            else:
                flag = False
        else:
            if x_index % 2 == 0:
                flag = False
            else:
                flag = True

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
            x_cart (float): the cartesian x coordinate
            y_cart (float): the cartesian y coordinate
            z_cart (float): the cartesian z coordinate
        """
        x_cart = (x_frac * x_res) + x_trans
        y_cart = (y_frac * y_res) + y_trans
        z_cart = (z_frac * z_res) + z_trans

        return x_cart, y_cart, z_cart

    return fractional_to_cartesian_coordinates

def create_cube_coords(x_index, y_index, z_index, frac_to_cart, coords):
    '''
    Caluculates the next 8 coords for the next volxel that has been thresholded
    previously (logic in loop that calls this one).
    Args:
        x_index (int)           Indecies to the x coord in the coords array.
        y_index (int)           Indecies to the y coord in the coords array.
        z_index (int)           Indecies to the z coord in the coords array.
        frac_to_cart (float*3=>float*3): fractional to cartesian conversin function
        coords (float list): The coordinates of the thresholded voxels passed by
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





def write_tets_to_files(nvoxel, points, cells, tet_array, output_file, ffea_out, vtk_out):
    """
    outupt the files
    Args:
        nvoxel (int): voxel count
        points (numpy.ndarray): coordinates of vertices (no duplicates)
        cells (numpy.ndarray): nvoxel by 8 array listing indices of
                               vertices in points for each voxel
        tet_array (int np.ndarray): 2D number of tets by four, the entry for each tet
                                    is a list of its four vertices in the points array
        output_file (pothlib.Path): name stem for ouput files
        ffea_out (bool): if true write ffea input files
        vtk_out (bool): if true write a vtk file
    Returns:
        None
    """
    # make the vtk tet connectivity
    cells_con = make_vtk_cell_connectivity(tet_array, len(cells))

    # make the gris (vtk scene)
    vtk_pts = vtk.vtkPoints()
    vtk_pts.SetData(vtk.util.numpy_support.numpy_to_vtk(points, deep=True))
    grid = vtk.vtkUnstructuredGrid() #create unstructured grid
    grid.SetPoints(vtk_pts) #assign points to grid
    grid.SetCells(vtk.VTK_TETRA, cells_con) #assign tet cells to grid

    #write vtk file
    if vtk_out:
        vtk_output(grid, output_file)

    # write tetgen file for ffea input
    if ffea_out:
        ffea_output(grid, points, output_file, nvoxel, tet_array)

def vtk_output(grid, output_file):
    """
    setup and use vtk writer
    Args:
        grid (vtk.vtkUnstructuredGrid): vtk scene
        output_file (pathlib.Path)
    """
    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetFileName(str(output_file.with_suffix(".vtk")))
    writer.SetInputData(grid)
    writer.Update()
    writer.Write()

def get_vtk_surface(grid):
    """
    get the surface polygons from a vtk scene.
    Args:
        grid (vtk.vtkUnstructuredGrid): vtk scene
    Return
        (float np.array): the vertices on the surface
        (vtkmodules.vtkCommonCore.vtkIdTypeArray) connectivity of surface polygons
        (int) numer of cells in surface
    """
    # make a surface filter to extract the geometric boundary
    surf_filt = vtk.vtkDataSetSurfaceFilter()
    surf_filt.SetInputData(grid)
    surf_filt.Update()

    # get the geometric boundary (the suface of the volume)
    surf = surf_filt.GetOutput()

    # get the points in the surface geomatry
    surf_points=np.array(surf.GetPoints().GetData())

    # get the surface polygons
    cells = surf.GetPolys()
    cell_count = cells.GetNumberOfCells()
    cell_data = cells.GetData()

    return surf_points, cell_data, cell_count

def ffea_output(grid, points, output_file, nvoxel, tet_array):
    """
    construct the faces and output the ffea input files
    Args:
        grid (vtk.vtkUnstructuredGrid): vtk scene
        points (float np.ndarray): duplicate free list of vertices
        output_file (pathlib.Path): name stem of output files
        nvoxel (int): the number of voxels
        tet_array (int np.ndarray): 2D number of tets by four, the entry for each tet
                                    is a list of its four vertices in the points array
    Returns:
        None
    """
    surf_points, cell_data, cell_count = get_vtk_surface(grid)

    # make a connectivity into the tets points
    original_ids = np.zeros((len(surf_points),), dtype='int16')
    for pos, point in enumerate(surf_points):
        # index of sourface point in the tet's points array
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
                      nvoxel,
                      tet_array,
                      points,
                      faces,
                      original_ids,
                      f'# created by {getpass.getuser()} on {date}')



def make_voxel_connectivity(nvoxel, coords):
    """
    construct a duplicate free list of vertices coordinates and a
    connectivity list mapping voxles to lists of eight vertices
    Args:
        nvoxel (int):                  the number of voxels
        coords (float numpy.ndarray):  an 2d array with each element being an array of length 8
                                       (one element for each vertex of a voxel and is the index
                                        of a value in the coords array) and the length
                                        of the overall array is the number of voxels that have
                                        been thresholded.
    Returns:
        points (float numpy.ndarray):  duplicate free list of voxel vertices
        cells (int numpy.ndarray):     nvoxel by 8 array listing indices of
                                       vertices in points for each voxel

    """
    # make unique list of vertices for indexing purposes
    points = np.unique(coords, axis=0)

    # create empty connectivity array, of vertex indices into the points list
    connectivity = np.zeros((nvoxel*8,), dtype='int16')

    # for each vertex in the orginal array the connectivity of
    # that index is assigned to the index of that vertex in the points
    for vertex_index, _ in enumerate(coords):
        point = coords[vertex_index]
        connectivity[vertex_index] = np.where((points==point).all(axis=1))[0]

    # convert connectivity to array of length nvoxel in which each entry
    # is an array of eight indices into the points array; the eight points
    # represent the eight corners of the voxel
    cells = np.resize(connectivity, (nvoxel,8))

    return points, cells

def make_tet_connectivity(nvoxel, cells, alternate):
    """
    convert voxel connectivity to tetrohedron connectivity
    Args:
        nvoxel (int): number of voxels
        cells (int numpy.ndarray): nvoxel by 8 array listing indices of
                                    vertices in points for each voxel
        alternate (int numpy.ndarray): indexed by voxel number, zero if voxel is
                                        odd else voxel is even
    Returns:
        tet_array (int np.ndarray): 2D number of tets by four, the entry for each tet
                                    is a list of its four vertices in the points array
    """
    tet_array = np.zeros((nvoxel*5,4), dtype='int16') #tet array for .ele

    #iterate over cubes and convert to tets
    for i, cube in enumerate(cells):
        if alternate[i]:
            connectivity_one_vox = odd_cube_tets(cube)
        else:
            connectivity_one_vox = even_cube_tets(cube)

        for tet_index, con_one_tet in enumerate(connectivity_one_vox):
            tet_array[(i*5)+tet_index] = con_one_tet

    return tet_array

def make_vtk_tet_connectivity(connectivities):
    """
    setup and use vtk writer
    Args:
        connectivities (int np.ndarray): 2D number of tets by four, the entry for each tet
                                         is a list of its four vertices in the points array
    Returns:
        (vtk.vtkCellArray): array holding the tet's connectivity as vtk data
    """
    #create vtk array for holding cells
    cells_con = vtk.vtkCellArray()

    # add tetrahedron cells to array
    for tet in connectivities:
        tetra = vtk.vtkTetra()
        for i, coord in enumerate(tet):
            tetra.GetPointIds().SetId(i, coord)

        #add tet data to vtk cell array
        cells_con.InsertNextCell(tetra)

    return cells_con

def make_vtk_cell_connectivity(tet_array, cell_count):
    """
    setup and use vtk writer
    Args:
        tet_array (int np.ndarray): 2D number of tets by four, the entry for each tet
                                    is a list of its four vertices in the points array
        cell_count (int): the number of voxels
    Returns:
        (vtk.vtkCellArray): array holding the tet's connectivity as vtk data
    """
    #create vtk array for holding cells
    cells_con = vtk.vtkCellArray()

    # add tetrahedron cells to array
    for i in range(cell_count):
        for tet in range(5):
            tet_con = tet_array[(i*5) + tet]

            tetra = vtk.vtkTetra()
            tetra.GetPointIds().SetId(0, tet_con[0])
            tetra.GetPointIds().SetId(1, tet_con[1])
            tetra.GetPointIds().SetId(2, tet_con[2])
            tetra.GetPointIds().SetId(3, tet_con[3])

            cells_con.InsertNextCell(tetra) #add tet data to vtk cell array

    return cells_con
