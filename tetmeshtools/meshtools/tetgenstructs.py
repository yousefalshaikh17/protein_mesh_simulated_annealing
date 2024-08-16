"""
Classes for data from a TetGen output file
format (https://wias-berlin.de/software/index.jsp?id=TetGen&lang=1).

--------------------------------

Licensed under the GNU General Public License, Version 3.0 (the "License"); you
may not use this file except in compliance with the License. You may obtain a
copy of the License at <https://www.gnu.org/licenses/gpl-3.0.html>.

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

This work was funded by Joanna Leng's EPSRC funded RSE Fellowship (EP/R025819/1)

@copyright 2023
@author: j.h.pickering@leeds.ac.uk and j.leng@leeds.ac.uk
"""
# set up linting conditions
# pylint: disable = import-error
import collections
import numpy as np

from tetmeshtools.vector3 import Vector3

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
