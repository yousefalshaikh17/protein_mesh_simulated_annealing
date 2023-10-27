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

 @author: jonathan pickering, joanna leng 08 Feb 23
"""
import numpy as np

class Vector3(list):
    """
    a 3d vector class extending list and providing
    doc, cross, normalize and is_close methods
    """

    @classmethod
    def make_from_list(cls, source):
        """
        constructor from list
        Args:
            source ([float]): list of at least 3 floats
        """
        return cls(source[0], source[1], source[2])

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
        return Vector3.make_from_list(np.cross(self, rhs))

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
        return Vector3.make_from_list(tmp)

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
