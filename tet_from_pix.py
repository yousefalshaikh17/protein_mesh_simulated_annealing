#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 11:04:41 2021

@author: mollygravett
"""

#print(mesh.cells[0])

import numpy as np
import vtk
import math
from vtk.util import numpy_support

file = 'clipped_fewcubes_map_thresholded_0p00202.vtk'

reader = vtk.vtkDataSetReader()
reader.SetFileName(file)
reader.Update()
vtk_data =  reader.GetOutput()
points = np.array(reader.GetOutput().GetPoints().GetData())

cells = np.array(reader.GetOutput().GetCells().GetData())
print(np.size(cells))
#cells_ = np.array([[4, 0, 1, 2, 3], [4, 0, 1, 2, 3]])

x = int((np.size(cells))/5)
new_cells=np.resize(cells,(x,5))

#print(new_cells[0][4])

# longest_edge_array=np.zeros((x,))
# surface_area_array=np.zeros((x,))

# for i in range(np.shape(new_cells)[0]):
#     tet=new_cells[i]
    
#     #locate points
#     point0 = points[tet[1]]
#     point1 = points[tet[2]]
#     point2 = points[tet[3]]
#     point3 = points[tet[4]]
    
#     #calculate distances
#     distance_p0_p1 = math.sqrt(np.square(point0-point1).sum(axis=0))
#     distance_p1_p2 = math.sqrt(np.square(point1-point2).sum(axis=0))
#     distance_p2_p0 = math.sqrt(np.square(point2-point0).sum(axis=0))
#     distance_p2_p3 = math.sqrt(np.square(point2-point3).sum(axis=0))
#     distance_p3_p1 = math.sqrt(np.square(point3-point1).sum(axis=0))
#     distance_p3_p0 = math.sqrt(np.square(point3-point0).sum(axis=0))
    
#     #triangles sorted in descending order of side length
#     base = -np.sort(-np.array([distance_p0_p1, distance_p1_p2, distance_p2_p0]))
#     right = -np.sort(-np.array([distance_p0_p1, distance_p3_p1, distance_p3_p0]))
#     left = -np.sort(-np.array([distance_p3_p0, distance_p2_p3, distance_p2_p0]))
#     back = -np.sort(-np.array([distance_p3_p1, distance_p1_p2, distance_p2_p3]))
#     longest_edge=-np.sort(-np.concatenate((base, right, left, back)))[0] #longest edge in tet
#     longest_edge_array[i]=longest_edge
    
#     #surface area
#     base_area = 0
#     right_area = 0
#     left_area = 0
#     back_area = 0
#     #dictionary = {'base_': base, ''}
#     for triangle in [base,right,left,back]:
#         area = (1/4)*math.sqrt((triangle[0] + (triangle[1] + triangle[2]))*(triangle[2] - (triangle[0] - triangle[1]))*(triangle[2] + (triangle[0] - triangle[1]))*(triangle[0] + (triangle[1] - triangle[2])))
#        # dictionary(triangle)+'area'.append()
#         if triangle.all() == base.all():
#             base_area=base_area+area
#         if triangle.all() == right.all():
#             right_area=right_area+area
#         if triangle.all() == left.all():
#             left_area=left_area+area
#         if triangle.all() == back.all():
#             back_area=back_area+area
#     surface_area = base_area + right_area + left_area + back_area
#     surface_area_array[i]=surface_area
# #print(surface_area_array)

# #writing to vtk
# surface_area_array=numpy_support.numpy_to_vtk(surface_area_array)
# surface_area_array.SetName("SURFACE_AREA")
# vtk_data.GetCellData().AddArray(surface_area_array)
# longest_edge_array=numpy_support.numpy_to_vtk(longest_edge_array)
# longest_edge_array.SetName("LONGEST_EDGE")
# vtk_data.GetCellData().AddArray(longest_edge_array)
# writer = vtk.vtkUnstructuredGridWriter()
# writer.SetFileName("Output.vtk")
# writer.SetInputData(vtk_data)
# writer.Update()
# writer.Write()


    
