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
from ffeamesh import coord_utility
from ffeamesh.ffea_write import write_ffea_output


def make_vtk_tet_connectivity(connectivities):
    """
    Setup tet connectivity for tetrahedron cells for the vtk writer. Duplicates
    are not removed before calling the function.
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
    Setup cell connectivity for hexahedron (cubes) cells for the vtk writer.
    Removes duplicate coordinates.
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

def get_vtk_surface(grid):
    """
    Get the surface polygons from a vtk scene.
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
