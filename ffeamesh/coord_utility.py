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

    Authors: Joanna Leng, Jonathan Pickering - University of Leeds
    Emails: J.Leng@leeds.ac.uk, J.H.Pickering@leeds.ac.uk
"""
# set up linting
# pylint: disable = import-error


from itertools import count
from collections import (namedtuple, OrderedDict)
import numpy as np



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