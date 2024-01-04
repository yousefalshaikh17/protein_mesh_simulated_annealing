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
import unittest
import tetmeshtools.tetprops as tp
from tetmeshtools.tetmeshtools.tetgenstructs import NodePoint

class TestTetProps(unittest.TestCase):
    """
    tests of Tetgen file reader
    """

    def setUp(self):
        """
        build a full test class
        """
        self.verts = [[0.0, 0.0, 0.0],
                      [0.0, 1.0, 0.0],
                      [0.0, 0.0, 1.0],
                      [1.0, 0.0, 0.0]]

        self.nodes = [NodePoint(0, x[0], x[1], x[2]) for x in self.verts]

    def tearDown(self):
        """
        clean up
        """
        pass

    def test_tet_volume(self):
        """
        test the tetrahedron volume calculation
        """
        vol = tp.tet_volume(self.nodes)
        self.assertAlmostEqual(vol, 0.1666666666, msg="tetprops.tet_volume")

    def test_tet_area(self):
        """
        test the tetrahedron area calculation
        """
        area = round(tp.tet_area(self.nodes), 3)
        self.assertAlmostEqual(area, 2.366, msg="tetprops.tet_area")

    def test_triangle_area(self):
        """
        test the triangle area calculation
        """
        area = tp.area_of_triangle(self.nodes)
        self.assertAlmostEqual(area, 0.50000, msg="tetprops.triangle_area")

    def test_edges_to_area_ratio_squared(self):
        """
        test the calculation of perimiter length to area
        """
        ratio = tp.edges_to_area_ratio_squared(self.nodes)
        self.assertAlmostEqual(ratio, 23.31370849, msg="tetprops.edges_to_area_ratio_squared")
