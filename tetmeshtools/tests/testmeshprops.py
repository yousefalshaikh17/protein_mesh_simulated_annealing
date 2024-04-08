"""
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
# set up linting
# pylint: disable = import-error

import unittest
import numpy as np
import tetmeshtools.tetprops as tp
import tetmeshtools.meshtools.tetgenstructs as tgs
import tetmeshtools.vector3 as v3

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
