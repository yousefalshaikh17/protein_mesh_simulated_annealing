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
