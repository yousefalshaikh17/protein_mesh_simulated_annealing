"""
a three dimensional vector

----------------------------------

Licensed under the GNU General Public License, Version 3.0 (the "License"); you
may not use this file except in compliance with the License. You may obtain a
copy of the License at <https://www.gnu.org/licenses/gpl-3.0.html>.

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

This work was funded by Joanna Leng's EPSRC funded RSE Fellowship (EP/R025819/1)

@copyright 2020
@author: j.h.pickering@leeds.ac.uk and j.leng@leeds.ac.uk
"""
import numpy as np

class Vector3(list):
    """
    a 3d vector class extending list and providing
    doc, cross, normalize and is_close methods
    """

    @classmethod
    def make_from_np_array(cls, source):
        """
        constructor from list
        Args:
            cls (<class>): class specifier
            source ([float]): list of at least 3 floats
        Returns
            Vector3
        """
        return cls([source[0], source[1], source[2]])

    def dot(self, rhs):
        """
        dot (inner) product self.rhs
        Args:
            rhs (NodePoint): right hand side of dot
        Returns
            float
        """
        return np.dot(self, rhs)

    def cross(self, rhs):
        """
        cross product self x rhs
        Args:
            rhs (NodePoint): right hand side of cross
        Returns:
            (Vector3)
        """
        return Vector3.make_from_np_array(np.cross(self, rhs))

    def length(self):
        """
        find length of vector
        Returns:
            (float): the lenght
        """
        return np.linalg.norm(self)

    def normalize(self):
        """
        make a normalized copy of objects vector
        Returns
            [float]
        """
        tmp = np.divide(self, self.length())
        return Vector3.make_from_np_array(tmp)

    def isclose(self, rhs):
        """
        floating point comparison of two vectors
        Args:
            rhs (Vector3): comaprison object
        Returns:
            True if all comonants pass numpy isclose
        """
        results = np.isclose(self, rhs)
        return all(x for x in results)
