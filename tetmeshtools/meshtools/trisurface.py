"""
container for a surface made of triangles

------------------------------------

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
# pylint: disable = import-error
import numpy as np

from tetmeshtools.meshtools.linesegment import LineSegment
import tetmeshtools.tetprops as tp

class TriSurface():
    """
    container for a surface made of triangles, defined by int reference to elements array
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
        # make list in which each node appears only once
        index_set = set()
        for index in self._nodes:
            index_set.add(index)

        surface_nodes = {}
        for index in index_set:
            surface_nodes[index] = self._nodes[index]

        return surface_nodes

    def get_surface_node_indices(self):
        """
        get the indices of the nodes that are used in the surface
        Returns:
            [int]
        """
        index_set = set()
        for index in self._nodes:
            index_set.add(index)

        return list(index_set)

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

    def get_triangle_signed_area_list(self):
        """
        get a list of the signed areas of the triangles
        """
        areas = []

        for index in self._faces:
            face_nodes = self.get_triangle_nodes(index)
            areas.append(tp.area_of_triangle(face_nodes, True))

        return areas

    def get_model_limits(self):
        """
        get the min and max for the axis
        Returns
            (<min x>, <max x>, <min y>, <max y>, <min z>, <max z>)
        """
        nodes = self.get_surface_nodes().values()

        x_values = [node.x for node in nodes]
        min_x = min(x_values)
        max_x = max(x_values)

        y_values = [node.y for node in nodes]
        min_y = min(y_values)
        max_y = max(y_values)

        z_values = [node.z for node in nodes]
        min_z = min(z_values)
        max_z = max(z_values)

        return min_x, max_x, min_y, max_y, min_z, max_z

    def __str__(self):
        """to string"""
        return f"Surface: {len(self._faces)} triangles on {len(self._nodes)} nodes"
