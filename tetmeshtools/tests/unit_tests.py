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

from tetmeshtools.tests.testreadtetgen import TestReadTetgen
from tetmeshtools.tests.testmeshprops import TestMeshProps
from tetmeshtools.tests.testtetprops import TestTetProps

def make_suite():
    """
    make a unittest TestSuite object
        Returns
            (unittest.TestSuite)
    """
    suite = unittest.TestSuite()

    suite.addTest(TestReadTetgen('test_read_nodes'))
    suite.addTest(TestReadTetgen('test_read_tets'))
    suite.addTest(TestReadTetgen('test_read_faces'))
    suite.addTest(TestMeshProps('test_tet_volume'))
    suite.addTest(TestMeshProps('test_nodepoint_to_edge_array'))
    suite.addTest(TestMeshProps('test_vector3'))
    suite.addTest(TestTetProps('test_tet_volume'))
    suite.addTest(TestTetProps('test_tet_area'))
    suite.addTest(TestTetProps('test_triangle_area'))
    suite.addTest(TestTetProps('test_edges_to_area_ratio_squared'))

    return suite

def run_all_tests():
    """
    run all tests in the TestSuite
    """
    runner = unittest.TextTestRunner()
    runner.run(make_suite())

if __name__ == '__main__':
    run_all_tests()
