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
