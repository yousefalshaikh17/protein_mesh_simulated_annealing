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
import numpy as np

import ffeamesh.tetmeshtools.tetgenread as tr

class TetMesh():
    """
    mutable store of a tetmesh derived from tetgen files
    """

    def __init__(self, nodes, tets):
        """
        set up the object
        Args:
            nodes (tetview.NodePoint)
            tets  (tetview.Tetrahedron4)
        """
        ## the vertices of the tets
        self._nodes = nodes

        ## the connectivity
        self._tets  = tets

    def get_nodes(self):
        """getter for the nodes"""
        return self._nodes

    def number_of_tets(self):
        """getter for the number of tets in mesh"""
        return len(self._tets)

    def get_tet_keys(self):
        """getter for the tet's dict keys (tetgen indices)"""
        return self._tets.keys()

    def get_tets(self):
        """getter for the tets"""
        return self._tets

    def get_tet_verts(self, index):
        """
        get the 4 nodes forming a tet
        Args:
            index (int): tetgen index
        Returns:
            [NodePoint]
        """
        tet = self._tets[index]

        verts = []
        verts.append(self._nodes[tet.vert0])
        verts.append(self._nodes[tet.vert1])
        verts.append(self._nodes[tet.vert2])
        verts.append(self._nodes[tet.vert3])

        return verts

    def get_tet_ctr(self, index):
        """
        get the center of the tetrahedron
        Args:
            index (int): array index of tet
        Returns:
            [float]: coordinates of center
        """
        total = [0.0, 0.0, 0.0]

        for vert in self.get_tet_verts(index):
            total = np.add(total, vert)

        return [x/4 for x in total]

    def get_tet_signed_6volume(self, index):
        """
        compute a six times the signed volume for a tet using the formular
        a, b, c are three edge vectors extending from one vertex
        V = ((axb).c)

        to convert to actual volume: vol = np.abs(v/6)
        """
        nodes = self.get_tet_verts(index)
        v01 = self._nodes[nodes[0].index].to_edge_array(self._nodes[nodes[1].index])
        v02 = self._nodes[nodes[0].index].to_edge_array(self._nodes[nodes[2].index])
        v03 = self._nodes[nodes[0].index].to_edge_array(self._nodes[nodes[3].index])

        tmp = np.cross(v01, v02)
        tmp = np.dot(tmp, v03)

        return tmp

    def get_tet_signed_6volume_list(self):
        """
        get the six times the signed volumes of the all the tets
        to convert to actual volume: vol = [np.abs(v/6) for v in vols6]
        Returns:
            [float]
        """
        sv = []
        for index in self._tets.keys():
            sv.append(self.get_tet_signed_6volume(index))

        return sv

    def __str__(self):
        """to string"""
        return f"Mesh: {len(self._tets)} tets on {len(self._nodes)} nodes"
