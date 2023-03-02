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
