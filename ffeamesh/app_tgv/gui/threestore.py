"""
Created on 22 Dec 2022

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this file except in compliance with the License. You may obtain a copy of the
License at http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

This work was funded by Joanna Leng's EPSRC funded RSE Fellowship (EP/R025819/1)

@copyright 2022
@author: j.h.pickering@leeds.ac.uk and j.leng@leeds.ac.uk
"""
import numbers
import numpy as np

class ThreeStore():
    """
    storage for three properties with axis name useful in graphics
    """
    def __init__(self, x=None, y=None, z=None):
        """
        initialize
        Args:
            x (any): value related to x axis
            y (any): value related to y axis
            z (any): value related to z axis
        """
        ## property for x
        self.x = x

        ## property for y
        self.y = y

        ## property for z
        self.z = z

    def length(self):
        """
        if the data is numeric return the 'length'
        """
        if isinstance(self.x, numbers.Number):
            if isinstance(self.y, numbers.Number) and isinstance(self.z, numbers.Number):
                return np.sqrt(self.x*self.x + self.y*self.y + self.z*self.z)

        return None

    def rotate(self, rot_x, rot_y):
        """
        rotate about
        """
        radius = self.length()

        if radius is None:
            raise TypeError("Store is not numeric and cannot be rotated")

        if np.isclose(0.0, radius):
            return self

        phi = np.radians(rot_x)
        theta = np.radians(rot_y)

        new_x = radius * np.sin(phi) * np.cos(theta)
        new_y = radius * np.sin(phi) * np.sin(theta)
        new_z = radius * np.cos(phi)

        return ThreeStore(new_x, new_y, new_z)

    def copy(self, other):
        """
        copy and object that has .x .y and .z
        """
        self.x = other.x
        self.y = other.y
        self.z = other.z

    def reverse(self):
        """
        reverse signs
        """
        return ThreeStore(-self.x, -self.y, -self.z)

    def __sub__(self, rhs):
        """
        provide subtraction
        """
        return ThreeStore(self.x - rhs.x,
                          self.y - rhs.y,
                          self.z - rhs.z)

    def __str__(self) -> str:
        """
        provide to string function
        """
        return f"<x:{self.x}, y:{self.y}, z:{self.z}>"

    def to_list(self):
        """
        return list of [x, y, z]
        """
        return [self.x, self.y, self.z]
