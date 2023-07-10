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
import collections
import numpy as np

from ffeamesh.vector3 import Vector3

## data structure for metadata line of tetgen .node file
## number of points, dimension, number of attributes and boundary markers
NodeMetaData = collections.namedtuple("NodeMetaData", "points, dimension, attributes, bms")

## 3D point with index base
_NodePoint = collections.namedtuple("_NodePoint", "index, x, y, z")

class NodePoint(_NodePoint):
    """
    extend 3D point
    """

    def to_array(self):
        """
        convert to numpy.array
        Return:
            numpy.array(float32)
        """
        return np.array([self.x, self.y, self.z], dtype=np.float32)

    def to_edge_array(self, rhs):
        """
        vector from self to rhs
        Args:
            rhs (NodePoint): the 'to' vector
        Returns:
            list(float)
        """
        return [rhs.x - self.x, rhs.y - self.y, rhs.z - self.z]

    def to_edge_vector(self, rhs):
        """
        vector from self to rhs
        Args:
            rhs (NodePoint): the 'to' vector
        Returns:
            EdgeVector
        """
        return Vector3(self.to_edge_array(rhs))

## data structure for metadata line of tetgen .node file
## number of faces and boundary markers
FaceMetaData = collections.namedtuple("FaceMetaData", "faces, bm")

## a face
Face = collections.namedtuple("Face", "index, vert0, vert1, vert2, bm")

## data structure for metadata line of tetgen .node file
## number of tets, number of nodes pre tet and region attribute
TetMetaData = collections.namedtuple("TetMetaData", "tets, nodes, ra")

## a tetrahedron
Tetrahedron4 = collections.namedtuple("Tetrahedron4", "index, vert0, vert1, vert2, vert3, ra")
