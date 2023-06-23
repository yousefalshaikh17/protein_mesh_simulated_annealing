#!/usr/bin/env python
"""
 six_tets_demo.py
 
 Create a pair of voxel with either 5 or 6 tetrahedron per
 voxel.
 
 
---------------------------------

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
import pathlib
import argparse
import numpy as np
import vtk.util.numpy_support

import ffeamesh.vtk_utility as vtk_u
from ffeamesh import vtk_write

"""
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

def unit_cube_coords():
    """vertices of unit cube"""
    return [[0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0], # cube 0 end
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 0.0, 1.0],
            [2.0, 0.0, 1.0],
            [2.0, 1.0, 1.0],
            [1.0, 1.0, 1.0]]

def cube_5tet_indices():
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

def cube_6tet_indices():
    """
    return a list of lists for the constuction of 5 tets from the 8 vertices of an even cube
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
            [0, 6, 4, 5],
            [8, 14, 13, 9],
            [8, 14, 9, 10],
            [8, 14, 10, 11],
            [8, 14, 11, 15],
            [8, 14, 15, 12],
            [8, 14, 12, 13]]

def get_args():
    """
    get the command line options
    """
    parser = argparse.ArgumentParser("make a pair of voxels")

    parser.add_argument("-5",
                        "--five_tets",
                        action='store_true',
                        help="if set decomp into 5 tets")

    return parser.parse_args()

def main():
    """run the script"""
    args = get_args()

    connectivity = None
    file_name = "cube"
    if args.five_tets:
        connectivity = cube_5tet_indices()
        file_name += f"_{5}tets.vkt"
    else:
        connectivity = cube_6tet_indices()
        file_name += f"_{6}tets.vkt"

    points_np = np.array(unit_cube_coords())
    cells_con = vtk_u.make_vtk_tet_connectivity(connectivity)

    # make the grid (vtk scene)
    vtk_pts = vtk.vtkPoints()
    vtk_pts.SetData(vtk.util.numpy_support.numpy_to_vtk(points_np, deep=True))
    grid = vtk.vtkUnstructuredGrid() #create unstructured grid
    grid.SetPoints(vtk_pts) #assign points to grid
    grid.SetCells(vtk.VTK_TETRA, cells_con) #assign tet cells to grid

    # write vtk file
    vtk_write.vtk_output(grid, pathlib.Path(file_name))

if __name__ == '__main__':
    main()
