#!/usr/bin/env python
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

Generating tetrahedral meshes from pixel data
"""
#from ffeamesh.mrc_zoom import mrc_zoom
import sys
import datetime
import argparse
import pathlib
import numpy as np
#import vtk
import mrcfile
import vtk.util.numpy_support
from ffeamesh.writers import write_ffea_output

#import vtkwmtk from vmtk
# from chimerax.map_data import mrc
# ## #coarsen
# g = mrc.open('map_equilibrium_mystructrigor_15A_0p00202.mrc')[0]


#to coarsen a mesh (to 15 Å for example) in chimerax use:
#vol resample #1 spacing 15 - this coarsens to 15 Å voxels
#save newmap.mrc model #2 - this saves your new coarsened model as "newmap.mrc"

def get_args():
    """
    get the command line arguments
    Returns:
        (argparse.namespace)
    """

    parser = argparse.ArgumentParser("""process MRC files and produces a regular
        tetrahedral volumetric mesh for FFEA using the "marching tet" algorithm.
        This is written out in the tetgen .ele, .face, and .node file format for
        later use in FFEA, and .vtk for mesh analysis.

        Coding:   Molly Gravett (bsmgr@leeds.ac.uk),
        Joanna Leng (J.Leng@leeds.ac.uk),
        Jarvellis Rogers (J.F.Rogers1@leeds.ac.uk)
        Jonathan Pickering (J.H.Pickering@leeds.ac.uk)""")

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

    parser.add_argument("-v",
                        "--vtk",
                        action="store_true",
                        help="produce vtk output.")

    parser.add_argument("-f",
                        "--ffea",
                        action="store_true",
                        help="produce ffea output.")

    parser.add_argument("-w",
                        "--overwrite",
                        action="store_true",
                        help="overwrite")

    return parser.parse_args()

def validate_command_line(args):
    """
    validate the command line arguments#
    Args:
        args (argeparse.Namespace): the command line arguments
    Returns
        (str): error message if fail else None
    """
    # check the input file can be found
    if not args.input.exists():
        return f"Error: file {args.input} does not exist!"

    # check for overwriteing files
    if not args.overwrite:
        if args.vtk:
            file = args.output.with_suffix(".vtk")
            if file.exists():
                return f"Error: file {file} exists, use option -w to allow overwrite."

        # TODO CHECK FOR FFEA FILES

    # check that some output has been specified
    if not args.vtk and not args.ffea:
        return "Error: you must specify and output type (vtk, ffea)"

    return None

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

def main():
    """
    run the script
    """
    args = get_args()
    error_message = validate_command_line(args)

    if error_message is not None:
        print(error_message, file=sys.stderr)
        return

    threshold = 0.0
    input_file = args.input
    output_file = args.output

    mrc = mrcfile.mmap(input_file, mode='r+')

    nvoxel = 0

    #find each point
    for z in range(0, mrc.header.nz):
        for y in range(0, mrc.header.ny):
            for x in range(0, mrc.header.nx):
                if mrc.data[z,y,x] > threshold:
                    nvoxel=nvoxel+1

    coords = np.zeros((nvoxel*8, 3))
    ncoord = 0
    res = float(mrc.header.cella['x'])/mrc.header.nx
    x_trans = float(mrc.header.origin['x'])
    y_trans = float(mrc.header.origin['y'])
    z_trans = float(mrc.header.origin['z'])
    alternate = np.zeros((nvoxel,))
    location = 0

    #create hex points
    for z in range(0, mrc.header.nz):
        for y in range(0, mrc.header.ny):
            for x in range(0, mrc.header.nx):
                if mrc.data[z,y,x] > threshold:
                    coord1 = [((x-0.5)*res)+x_trans, ((y-0.5)*res)+y_trans, ((z-0.5)*res)+z_trans]
                    coords[ncoord] = coord1
                    coord2 = [((x+0.5)*res)+x_trans, ((y-0.5)*res)+y_trans, ((z-0.5)*res)+z_trans]
                    coords[ncoord+1] = coord2
                    coord3 = [((x+0.5)*res)+x_trans, ((y+0.5)*res)+y_trans, ((z-0.5)*res)+z_trans]
                    coords[ncoord+2] = coord3
                    coord4 = [((x-0.5)*res)+x_trans, ((y+0.5)*res)+y_trans, ((z-0.5)*res)+z_trans]
                    coords[ncoord+3] = coord4
                    coord5 = [((x-0.5)*res)+x_trans, ((y-0.5)*res)+y_trans, ((z+0.5)*res)+z_trans]
                    coords[ncoord+4] = coord5
                    coord6 = [((x+0.5)*res)+x_trans, ((y-0.5)*res)+y_trans, ((z+0.5)*res)+z_trans]
                    coords[ncoord+5] = coord6
                    coord7 = [((x+0.5)*res)+x_trans, ((y+0.5)*res)+y_trans, ((z+0.5)*res)+z_trans]
                    coords[ncoord+6] = coord7
                    coord8 = [((x-0.5)*res)+x_trans, ((y+0.5)*res)+y_trans, ((z+0.5)*res)+z_trans]
                    coords[ncoord+7] = coord8
                    ncoord=ncoord+8
                    #logic for alternating tet division 0 (a) or 1 (b)
                    if z % 2 == 0:
                        if y % 2 == 0:
                            if x % 2 == 0:
                                alternate[location] = 0
                            else:
                                alternate[location] = 1
                        else:
                            if x % 2 == 0:
                                alternate[location] = 1
                            else:
                                alternate[location] = 0
                    else:
                        if y % 2 == 0:
                            if x % 2 == 0:
                                alternate[location] = 1
                            else:
                                alternate[location] = 0
                        else:
                            if x % 2 == 0:
                                alternate[location] = 0
                            else:
                                alternate[location] = 1
                    location=location+1

    points_ = np.unique(coords, axis=0) #make unique list of points for index

    #create connectivity of cubes using index from point list
    connectivity_ = np.zeros((nvoxel*8,), dtype='int16')
    for pos in range(len(coords)):
        point = coords[pos]
        connectivity_[pos] = np.where((points_==point).all(axis=1))[0]
    cells_ = np.resize(connectivity_, (nvoxel,8))

    cells_con = vtk.vtkCellArray() #create vtk cell array
    tet_array = np.zeros((nvoxel*5,4), dtype='int16') #tet array for .ele

    #iterate over cubes and convert to tets
    for i in range(len(cells_)):
        hex=cells_[i]
        cube = hex[0:8]
        if alternate[i] == 0:
            connectivity = even_cube_tets(cube)
        elif alternate[i] == 1:
            connectivity = odd_cube_tets(cube)
        for tet in range(len(connectivity)):
            tet_array[(i*5)+tet] = connectivity[tet]
            tet_con = connectivity[tet]
            tetra = vtk.vtkTetra()
            tetra.GetPointIds().SetId(0, tet_con[0])
            tetra.GetPointIds().SetId(1, tet_con[1])
            tetra.GetPointIds().SetId(2, tet_con[2])
            tetra.GetPointIds().SetId(3, tet_con[3])
            cells_con.InsertNextCell(tetra) #add tet data to vtk cell array


    vtkPts = vtk.vtkPoints()
    vtkPts.SetData(vtk.util.numpy_support.numpy_to_vtk(points_, deep=True))
    grid = vtk.vtkUnstructuredGrid() #create unstructured grid
    grid.SetPoints(vtkPts) #assign points to grid
    grid.SetCells(vtk.VTK_TETRA, cells_con) #assign tet cells to grid

    surfFilt2 = vtk.vtkDataSetSurfaceFilter()
    surfFilt2.SetInputData(grid)
    surfFilt2.Update()
    surf = surfFilt2.GetOutput()
    surf_points=np.array(surf.GetPoints().GetData())
    cells = surf.GetPolys()
    nCells = cells.GetNumberOfCells()
    array = cells.GetData()
    original_ids = np.zeros((len(surf_points),), dtype='int16')
    for pos in range(len(surf_points)):
        point = surf_points[pos]
        original_ids[pos] = np.where((points_==point).all(axis=1))[0]

    # This holds true if all polys are of the same kind, e.g. triangles.
    assert(array.GetNumberOfValues()%nCells==0)
    nCols = array.GetNumberOfValues()//nCells
    numpy_cells = np.array(array)
    faces = numpy_cells.reshape((-1,nCols))

    #write to vtk
    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetFileName(str(output_file.with_suffix(".vtk")))
    writer.SetInputData(grid)
    writer.Update()
    writer.Write()

    #write to tetgen .ele, .node, .face
    date = datetime.datetime.now().strftime("%x")
    comment = '# created by Molly Gravett ' + date
    write_ffea_output(output_file, nvoxel, tet_array, points_, faces, original_ids, comment)

if __name__ == "__main__":
    main()
