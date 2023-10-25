"""
@copyright 2022
@author: j.h.pickering@leeds.ac.uk and j.leng@leeds.ac.uk
"""
# set up linting conditions
# pylint: disable = import-error
import pathlib

import ffeamesh.tetmeshtools.tetgenwrite as tw

class TetModel():
    """
    class holding a complete tetgen mesh
    """
    def __init__(self, surface, tets, source=None):
        ## the surface object
        self._surface = surface

        ## the 3d model
        self._tets = tets

        ## data source
        self._source = source

    def get_source(self):
        """
        getter for the source
        """
        return self._source

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

    def get_tet_verts(self, index):
        """
        getter for the vertices of a tet
        Args:
            index (int)
        """
        return self._tets.get_tet_verts(index)

    def get_tet_nodes(self, index):
        """
        getter for the nodes forming a tet
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
        # TODO no method get_tet in class
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

    def get_tetmesh(self):
        """get the tetmesh object"""
        return self._tets

    def get_surface(self):
        """getter for the surface"""
        return self._surface

    def __str__(self):
        """to string method"""
        return f"TetModel: <{self._surface}> <{self._tets}>"

    def write_to_file(self, file_name):
        """
        write contents to file_name.1.ele, file_name.1.node, file_name.1.face
        Args:
            file_name (str)
        """
        path = pathlib.Path(file_name)

        tw.write_tetgen_nodes(path, self._surface.get_nodes())
        tw.write_tetgen_faces(path, self._surface.get_faces())
        tw.write_tetgen_elements(path, self._tets.get_tets())
