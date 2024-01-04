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
from collections import (namedtuple, OrderedDict)

## data structure for a 3D coordinate
Coordinate = namedtuple("Coordinate", "x, y, z")

class CoordTransform():
    """
    store for the fractional and cartesian coordinates of a 3D lattice point
    """
    def __init__(self, frac, cart):
        """
        initalize the object
        Args:
            frac (Coordinate): fractional coordinate
            cart (Coordinate): cartesian coordinate
        """
        ## storage for the fractional coordinate
        self.frac = frac

        ## storage for the cartesian coordinate
        self.cart = cart

    def __eq__(self, lhs):
        """
        equality is based on fractional coordinate only
        Args:
            lhs (CoordTransform): the object for comparison
        Returns:
            (bool) True if fractional coordinates are equal else False
        """
        if self.frac.x != lhs.frac.x:
            return False

        if self.frac.y != lhs.frac.y:
            return False

        if self.frac.z != lhs.frac.z:
            return False

        return True

    def __hash__(self):
        """
        hash value is based on the fractional coordinate only
        Returns:
            (int) the hash value of the fractional coordinate
        """
        return hash((self.frac.x, self.frac.y, self.frac.z))

    def __repr__(self):
        """
        string representation of the object
        """
        str_frac = f"({self.frac.x}, {self.frac.y}, {self.frac.z})"
        str_cart = f"({self.cart.x}, {self.cart.y}, {self.cart.z})"
        return f"<CoordTransform: {str_frac} => {str_cart}>"


class UniqueTransformStore():
    """
    store for CoordTransform which only holds unique objects,
    based on fractional coordinate
    """
    def __init__(self):
        """
        initalize the object
        """
        ## current size of the array
        self.current_size = 0

        ## dictionary holding the data, will work python < 3.7
        self.data = OrderedDict()

    def add(self, coord):
        """
        add a new CoordTransform
        Args:
            coord (CoordTransform): the object to be added
        Returns
            (int): the index of the object in the list of unique objects
        """
        # if already stored return the equivalent array index
        if coord in self.data.keys():
            return self.data[coord]

        # if new point store and return the equivalent array index
        self.data[coord] = self.current_size
        tmp = self.current_size
        self.current_size += 1

        return tmp

    def to_list(self):
        """
        convert the keys to a list in entry order
        Returns:
            (CoordTransform list)
        """
        return self.data.keys()

    def __str__(self):
        """
        provide string representation
        """
        return f"<UniqueArray: {self.current_size} items>"

    def __hash__(self):
        """
        hash function is the hash of the data
        """
        return hash(self.data)
