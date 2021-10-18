# -*- coding: utf-8 -*-
#
#  This file is part of the FFEA simulation package
#
#  Copyright (c) by the Theory and Development FFEA teams,
#  as they appear in the README.md file.
#
#  FFEA is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  FFEA is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
#
#  To help us fund FFEA development, we humbly ask that you cite
#  the research papers on the package.
#

"""
        tet_from_pix.py
        Authors: Molly Gravett, Joanna Leng, Jarvellis Rogers - University of Leeds
        Emails: bsmgr@leeds.ac.uk, J.Leng@leeds.ac.uk, J.F.Rogers1@leeds.ac.uk
"""

#Generating tetrahedral meshes from pixel data

#DEBUG - test command:
# python tet_from_pix.py --input="data/newmap_15A_0p00878.mrc" --output="new_0p00878_mesh_15A" --threshold=0.00878

from os import path
import numpy as np
import vtk
import mrcfile
import vtk.util.numpy_support
import datetime
import getopt
import sys
#import vtkwmtk from vmtk
# from chimerax.map_data import mrc
# ## #coarsen
# g = mrc.open('map_equilibrium_mystructrigor_15A_0p00202.mrc')[0]


#to coarsen a mesh (to 15 Å for example) in chimerax use:
#vol resample #1 spacing 15 - this coarsens to 15 Å voxels
#save newmap.mrc model #2 - this saves your new coarsened model as "newmap.mrc"

def pathCheck(fpath):
    if not path.exists(fpath):
        msg = "ERROR: Input file path does not exist. Please try again."
        sys.exit(msg)

try:
    options, remainder = getopt.getopt(sys.argv[1:], "i:o:t:", ["input=", "output=", "threshold="])
except getopt.GetoptError as err:
    print("ERROR: " + str(err) + "\n")
    usage()

for opt, arg in options:
    if opt in ("-i", "--input"):
        mrcfilename = arg
        pathCheck(mrcfilename)
    elif opt in ("-o", "--output"):
        chosen_filename = arg
    elif opt in ("-t", "--threshold"):
        try:
            threshold = float(arg)
        except ValueError:
            msg = "ERROR: Threshold value must be a number."
            sys.exit(msg)

# Checks to see if mandatory options have been called
try:
    mrcfilename
except NameError:
    msg = "ERROR: Input MRC file path not defined."
    sys.exit(msg)
try:
    chosen_filename
except NameError:
    msg = "ERROR: Output destination file path not defined."
    sys.exit(msg)
try:
    threshold
except NameError:
    msg = "ERROR: Threshold not defined."
    sys.exit(msg)

mrc = mrcfile.open(mrcfilename, mode='r+')

a = mrc.data

nx = mrc.header.nx
ny = mrc.header.ny
nz = mrc.header.nz

nvoxel = 0

#find each point
for z in range(0, nz):
    for y in range(0, ny):
        for x in range(0, nx):
            if a[z,y,x] >= threshold:
                nvoxel=nvoxel+1

coords = np.zeros((nvoxel*8, 3))
ncoord = 0
res = float(mrc.header.cella['x'])/nx
x_trans = float(mrc.header.origin['x'])
y_trans = float(mrc.header.origin['y'])
z_trans = float(mrc.header.origin['z'])
alternate = np.zeros((nvoxel,))
location = 0



#create hex points
for z in range(0, nz):
    for y in range(0, ny):
        for x in range(0, nx):
            if a[z,y,x] >= threshold:
                coord1 = [(x*res)+x_trans, (y*res)+y_trans, (z*res)+z_trans]
                coords[ncoord] = coord1
                coord2 = [((x+1)*res)+x_trans, (y*res)+y_trans, (z*res)+z_trans]
                coords[ncoord+1] = coord2
                coord3 = [((x+1)*res)+x_trans, ((y+1)*res)+y_trans, (z*res)+z_trans]
                coords[ncoord+2] = coord3
                coord4 = [(x*res)+x_trans, ((y+1)*res)+y_trans, (z*res)+z_trans]
                coords[ncoord+3] = coord4
                coord5 = [(x*res)+x_trans, (y*res)+y_trans, ((z+1)*res)+z_trans]
                coords[ncoord+4] = coord5
                coord6 = [((x+1)*res)+x_trans, (y*res)+y_trans, ((z+1)*res)+z_trans]
                coords[ncoord+5] = coord6
                coord7 = [((x+1)*res)+x_trans, ((y+1)*res)+y_trans, ((z+1)*res)+z_trans]
                coords[ncoord+6] = coord7
                coord8 = [(x*res)+x_trans, ((y+1)*res)+y_trans, ((z+1)*res)+z_trans]
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

#tet division for even cubes
def a(cube):
    tet1 = np.array([cube[0], cube[4], cube[5], cube[7]])
    tet2 = np.array([cube[0], cube[1], cube[2], cube[5]])
    tet3 = np.array([cube[2], cube[5], cube[6], cube[7]])
    tet4 = np.array([cube[0], cube[2], cube[3], cube[7]])
    tet5 = np.array([cube[0], cube[2], cube[5], cube[7]])
    tet_list = [tet1, tet2, tet3, tet4, tet5]
    return tet_list

#tet division for odd cubes
def b(cube):
    tet1 = np.array([cube[0], cube[1], cube[3], cube[4]])
    tet2 = np.array([cube[1], cube[4], cube[5], cube[6]])
    tet3 = np.array([cube[1], cube[2], cube[3], cube[6]])
    tet4 = np.array([cube[3], cube[4], cube[6], cube[7]])
    tet5 = np.array([cube[1], cube[3], cube[4], cube[6]])
    tet_list = [tet1, tet2, tet3, tet4, tet5]
    return tet_list

#iterate over cubes and convert to tets
for i in range(len(cells_)):
    hex=cells_[i]
    cube = hex[0:8]
    if alternate[i] == 0:
        connectivity = a(cube)
    elif alternate[i] == 1:
        connectivity = b(cube)
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
writer.SetFileName(chosen_filename+".vtk")
writer.SetInputData(grid)
writer.Update()
writer.Write()

#write to tetgen .ele, .node, .face
#.ele
#First line: <# of tetrahedra> <nodes per tetrahedron> <# of attributes>
#Remaining lines list of # of tetrahedra:
#<tetrahedron #> <node> <node> <node> <node> ... [attributes]
date = datetime.datetime.now().strftime("%x")
comment = '# created by Molly Gravett ' + date
ele = open(chosen_filename+".1.ele", "w")
ele_first = str(nvoxel*5)+' 4 0\n'
ele.write(ele_first)
for i in range(len(tet_array)):
    ele_next=str(i+1)+' '+str(tet_array[i][0]+1)+' '+str(tet_array[i][1]+1)+' '+str(+tet_array[i][2]+1)+' '+str(tet_array[i][3]+1)+'\n'
    ele.write(ele_next)
ele.write(comment)
ele.close()

#.node
#First line: <# of points> <dimension (must be 3)> <# of attributes> <# of boundary markers (0 or 1)>
#Remaining lines list # of points:
#<point #> <x> <y> <z>
node = open(chosen_filename+".1.node", "w")
node_first = str(len(points_))+' 3 0 0\n'
node.write(node_first)
for i in range(len(points_)):
    node_next=str(i+1)+' '+str(points_[i][0])+' '+str(points_[i][1])+' '+str(+points_[i][2])+'\n'
    node.write(node_next)
node.write(comment)
node.close()

#.face
#First line: <# of faces> <boundary marker (0 or 1)>
#Remaining lines list of # of faces:
#<face #> <node> <node> <node> [boundary marker]

face = open(chosen_filename+".1.face", "w")
face_first = str(len(faces))+' 1\n'
face.write(face_first)
for i in range(len(faces)):
    face_next=str(i+1)+' '+str(original_ids[faces[i][1]]+1)+' '+str(original_ids[faces[i][2]]+1)+' '+str(original_ids[faces[i][3]]+1)+' -1\n'
    face.write(face_next)
face.write(comment)
face.close()
