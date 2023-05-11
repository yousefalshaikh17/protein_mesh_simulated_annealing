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

@copyright 2023
@author: j.h.pickering@leeds.ac.uk and j.leng@leeds.ac.uk
"""
# pylint: disable = import-error
import numpy as np

from ffeamesh.tetmeshtools.linesegment import LineSegment
import ffeamesh.tetprops as tp

class TriSurface():
    """
    container for a surface made of triangles
    """

    def __init__(self, nodes, faces):
        """
        set up the object
        Args:
            nodes (NodePoint)
            faces (Face)
        """
        ## the vertices of the faces
        self._nodes = nodes

        ## the connectivity dict{index => tetgenstruct.Face}
        self._faces = faces

    def get_faces(self):
        """getter for the faces"""
        return self._faces

    def get_nodes(self):
        """getter for the nodes"""
        return self._nodes

    def get_edge_indices(self):
        """
        get array of face edges
        Return:
            [[int, int]]: start end pairs
        """
        pairs = set()
        for face in self._faces:
            start = face.vert0
            end = face.vert1
            pairs.add(LineSegment(start, end))
            start = end
            end = face.vert2
            pairs.add(LineSegment(start, end))
            start = end
            end = face.vert0
            pairs.add(LineSegment(start, end))

        return [x.to_list() for x in pairs]

    def get_surface_ogl(self):
        """
        get the surface as an array of homoginised 3D coords and normals
        [[vertex], [normal].....]
        """
        triangles = []
        for face in self._faces:
            node0 = self._nodes[face.vert0]
            node1  = self._nodes[face.vert1]
            node2  = self._nodes[face.vert2]

            edge01 = node0.to_edge_array(node1)
            edge12 = node1.to_edge_array(node2)
            normal = np.cross(edge01, edge12)
            normal = np.array([normal[0], normal[1], normal[2], 1.0], dtype=np.float32)

            triangles.append(node0.to_array(1.0))
            triangles.append(normal)
            triangles.append(node1.to_array(1.0))
            triangles.append(normal)
            triangles.append(node2.to_array(1.0))
            triangles.append(normal)

        return triangles

    def get_surface_ctr(self):
        """
        get the center of the surface
        Args:
            index (int): array index of tet
        Returns:
            [float]: coordinates of center
        """
        total = [0.0, 0.0, 0.0]
        surface_nodes = self.get_surface_nodes()

        for node in surface_nodes.values():
            total = np.add(total, node.to_array())

        return [x/len(self._nodes) for x in total]

    def get_surface_nodes(self):
        """
        getter for the nodes that appear in the surface
        Returns:
            dict index => node: the nodes used in the surface
        """
        index_set = set()
        for index in self._nodes:
            index_set.add(index)

        surface_nodes = {}
        for index in index_set:
            surface_nodes[index] = self._nodes[index]

        return surface_nodes

    def get_triangle_verts(self, index):
        """
        get the node indices for a triangle
        Args:
            index (int): index number of triangle
        Returns:
            [int, int, int]
        """
        return [self._faces[index].vert0,
                self._faces[index].vert1,
                self._faces[index].vert2]

    def get_triangle_nodes(self, index):
        """
        get the node for a triangle
        Args:
            index (int): index number of triangle
        Returns:
            [NodePoint, NodePoint, NodePoint]
        """
        face = self._faces[index]
        return [self._nodes[face.vert0],
                self._nodes[face.vert1],
                self._nodes[face.vert2]]

    def surface_area(self):
        """
        get the total area of the surface
        """
        total = 0.0

        for index in self._faces:
            face_nodes = self.get_triangle_nodes(index)
            total += tp.area_of_triangle(face_nodes)

        return total

    def __str__(self):
        """to string"""
        return f"Surface: {len(self._faces)} triangles on {len(self._nodes)} nodes"
