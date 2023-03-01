"""
@copyright 2022
@author: j.h.pickering@leeds.ac.uk and j.leng@leeds.ac.uk
"""
# pylint: disable = import-error
import numpy as np

from tetview.linesegment import LineSegment

class TriSurface():
    """
    container for a surface made of triangles
    """

    def __init__(self, nodes, faces):
        """
        set up the object
        Args:
            nodes (tetview.NodePoint)
            faces (tetview.Face)
        """
        ## the vertices of the faces
        self._nodes = nodes

        ## the connectivity
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
        getter for the nodes that appear in the
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
