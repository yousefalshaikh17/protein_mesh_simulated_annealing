"""
You should have received a copy of the GNU General Public License.
If not, see <http://www.gnu.org/licenses/>.

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

This work was funded by Joanna Leng's EPSRC funded RSE Fellowship (EP/R025819/1)

@copyright 2023
@author: j.h.pickering@leeds.ac.uk and j.leng@leeds.ac.uk
"""
import unittest
import tetmeshtools.tetprops as tp
from tetmeshtools.meshtools.tetgenstructs import NodePoint

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
