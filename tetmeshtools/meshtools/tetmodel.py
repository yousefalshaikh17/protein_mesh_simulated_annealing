"""
top level container for a tetgen mesh model, stores the both the
surface triangle mesh and the volume tetrahedron mesh

-----------------------------------

Licensed under the GNU General Public License, Version 3.0 (the "License"); you
may not use this file except in compliance with the License. You may obtain a
copy of the License at <https://www.gnu.org/licenses/gpl-3.0.html>.

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

This work was funded by Joanna Leng's EPSRC funded RSE Fellowship (EP/R025819/1)

@copyright 2022
@author: j.h.pickering@leeds.ac.uk and j.leng@leeds.ac.uk
"""
# set up linting conditions
# pylint: disable = import-error
import pathlib
import numpy as np

import tetmeshtools.meshtools.tetgenwrite as tw

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

    def get_ctr(self):
        """
        get the geometric ctr of the model
        """
        return self._surface.get_surface_ctr()

    def get_bounding_sphere(self):
        """
        get the bounding sphere of model
        Returns
            [float, float, float]: centre
            float: radius
        """
        ctr = self._surface.get_surface_ctr()

        tmp = self._surface.get_surface_nodes().values()
        lengths = [np.linalg.norm(np.subtract(node.to_array(), ctr)) for node in tmp]

        return ctr, max(lengths)

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

    def find_peak_nodes(self):
        """
        fill list with array indices of nodes beloning to only one tet
        """
        node_indices = self._tets.get_nodes().keys()
        inclusion_counts = {index:[] for index in node_indices}

        for index, tet in self._tets.get_tets().items():
            inclusion_counts[tet.vert0].append(index)
            inclusion_counts[tet.vert1].append(index)
            inclusion_counts[tet.vert2].append(index)
            inclusion_counts[tet.vert3].append(index)

        lengths = []
        for value in inclusion_counts.values():
            lengths.append(len(value))

        print(f" Min {min(lengths)}\n Max {max(lengths)}")

    def node_to_triang(self):
        """
        make map of nodes to faces in which they are included
        """
        # dict of tet indices to list of nodes used as vertices
        node_indices = self._surface.get_surface_nodes().keys()
        inclusion_counts = {index:[] for index in node_indices}

        for key, face in self._surface.get_faces().items():
            inclusion_counts[face.vert0].append(key)
            inclusion_counts[face.vert1].append(key)
            inclusion_counts[face.vert1].append(key)
