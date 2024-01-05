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
# TODO only used in cost function

class TriangleNodes():
    """
    storage for the node indices of a triangle
    """

    def __init__(self, vert0, vert1, vert2):
        """setup object"""
        self._verts = [vert0, vert1, vert2]

    def __hash__(self):
        """symmetric hash function (order independant)"""
        return hash(self._verts[0])^hash(self._verts[1])^hash(self._verts[2])

    def __contains__(self, item):
        """over ride in operator"""
        return item in self._verts

    def __eq__(self, rhs):
        """order independant equals"""
        for vert in self._verts:
            if vert not in rhs:
                return False

        return True

    def __repr__(self):
        """string description"""
        return f"TriangleNodes({self._verts[0]}, {self._verts[1]}, {self._verts[2]})"
