"""
@copyright 2022
@author: j.h.pickering@leeds.ac.uk and j.leng@leeds.ac.uk
"""
import numpy as np

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
