"""
@copyright 2022
@author: j.h.pickering@leeds.ac.uk and j.leng@leeds.ac.uk
"""

class TetModel():
    """
    class holding a complete tetgen mesh
    """
    def __init__(self, surface, tets):
        self._surface = surface
        self._tets = tets

    def get_surface_array3f(self):
        """
        get the surface array as list of lists of 3 float values
        """
        return None

    def get_model_ctr(self):
        """
        get the geometric ctr of the model
        """
        return self._surface.get_surface_ctr()

    def get_face_edges(self):
        """
        get array of face edges
        Return:
            [[float, float, float]]: start end pairs
        """
        return self._surface.get_face_edges()

    def get_vertices(self):
        """
        getter for the nodes
        """
        nodes = self._surface.get_nodes()

        largest = max(nodes.keys())+1
        vertices = [[0.0, 0.0, 0.0] for _ in range(largest)]

        for node in self._surface.get_nodes().values():
            vertices[node.index][0] = node.x
            vertices[node.index][1] = node.y
            vertices[node.index][2] = node.z

        return vertices

    def get_tet(self, index):
        """
        getter for the indices of a tet
        """
        return self._tets.get_tet_verts(index)

    def get_tet_display_list(self, index):
        """
        get the indices for displaying a tet
        Args:
            index (int): tet's index
        Returns
            [int]: points array
            [int]: edges array
        """
        points = []
        for point in self.get_tet(index):
            points.append(point.index)

        lines = []
        lines.append(points[0])
        lines.append(points[1])
        lines.append(points[0])
        lines.append(points[2])
        lines.append(points[0])
        lines.append(points[3])
        lines.append(points[1])
        lines.append(points[2])
        lines.append(points[2])
        lines.append(points[3])
        lines.append(points[3])
        lines.append(points[1])

        return points, lines

    def get_edges(self):
        """
        getter for the indices of the edges
        """
        return self._surface.get_edge_indices()

    def get_display_surface(self):
        """
        get the triangles and normals
        """
        return self._surface.get_surface_ogl()

    def surface_triangle_count(self):
        """
        getter for the number of triangles in the surface
        """
        return len(self._surface.get_faces())
