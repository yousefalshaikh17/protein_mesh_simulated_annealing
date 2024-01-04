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
class LineSegment():
    """
    data structure for an LineSegment defined in terms of node indices
    providing order dependant equality: (1, 2) not equal to (2, 1)
    """

    def __init__(self, start, end):
        """
        set-up
        Args:
            start (int): index of start node of segment
            end (int): index of end node of segment
        """

        ## start node
        self.start = start

        ## end node
        self.end = end

    def __eq__(self, rhs):
        """
        == oprerator
        Args:
            rhs (LineSegment): right hand side of operator
        """
        # strict equality
        if self.start == rhs.start and self.end == rhs.end:
            return True

        # reverse direction equality
        return self.start == rhs.end and self.end == rhs.start

    def __hash__(self) -> int:
        """
        hash function is hash of tuple so order dependant
        Returns:
            int: hash value
        """
        return hash((self.start, self.end))

    def __repr__(self):
        """
        string from which object could be reconstructed
        Returns:
            str: text form of object
        """
        return f"LineSegment(start={self.start}, end={self.end})"

    def __str__(self):
        """
        string describing object in conversational style
        Returns:
            str: text describing object
        """
        return f"line from {self.start} to {self.end}"

    def to_list(self):
        """
        return the object as a list
        Returns:
            [start, end]
        """
        return [self.start, self.end]
