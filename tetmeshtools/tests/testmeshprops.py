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

    Authors: Joanna Leng, Jonathan Pickering - University of Leeds
    Emails: J.Leng@leeds.ac.uk, J.H.Pickering@leeds.ac.uk
"""
# set up linting
# pylint: disable = import-error

import unittest
import numpy as np
import ffeamesh.tetprops as tp
import ffeamesh.tetmeshtools.tetgenstructs as tgs
import ffeamesh.vector3 as v3

class TestMeshProps(unittest.TestCase):
    """
    tests of Tetgen file reader
    """

    def setUp(self):
        """
        build a full test class
        """
        pass

    def tearDown(self):
        """
        clean up
        """
        pass

    def test_nodepoint_to_edge_array(self):
        """
        test NodePoint to_edge_arry
        """
        p1 = tgs.NodePoint(0, 1.0, 1.0, 1.0)
        p2 = tgs.NodePoint(0, 2.0, 2.0, 2.0)

        result = p1.to_edge_array(p2)
        self.assertAlmostEqual(result,
                               [1.0, 1.0, 1.0],
                               msg="tetmeshtools.tetgenstructs.to_edge_array")

    def test_tet_volume(self):
        """
        test calculating tet volumes
        (0,0,0) (18,0,0) (9,18,0) (9,9,9)
        """
        target = 486.0

        coords = []
        coords.append(tgs.NodePoint(0,  0.0,  0.0, 0.0))
        coords.append(tgs.NodePoint(0, 18.0,  0.0, 0.0))
        coords.append(tgs.NodePoint(0,  9.0, 18.0, 0.0))
        coords.append(tgs.NodePoint(0,  9.0,  9.0, 9.0))

        volume = tp.tet_volume(coords)
        self.assertEqual(volume, target, msg="tetprops.tet_volume")

    def test_vector3(self):
        """
        test vector ops
        """
        vec1 = v3.Vector3([1.0, 1.0, 1.0])
        vec1a = v3.Vector3([1.0, 1.0, 1.000000000000001])
        vec2 = v3.Vector3([1.0, -1.0, 0.0])
        root3 = np.sqrt(3.0)
        norm = v3.Vector3([1.0/root3]*3)

        self.assertAlmostEqual(vec1.dot(vec2), 0.0, msg='Vector3.dot')
        self.assertAlmostEqual(vec1.length(), root3, msg="Vector3.length")
        self.assertAlmostEqual(vec1.normalize(), norm, msg='Vector3.normalize')
        self.assertTrue(vec1.isclose(vec1a), msg='Vector3.is_close')
