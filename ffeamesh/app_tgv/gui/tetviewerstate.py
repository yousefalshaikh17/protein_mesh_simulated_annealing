"""
Created on 24 March 2023

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this file except in compliance with the License. You may obtain a copy of the
License at http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

This work was funded by Joanna Leng's EPSRC funded RSE Fellowship (EP/R025819/1)

@copyright 2023
@author: j.h.pickering@leeds.ac.uk and j.leng@leeds.ac.uk
"""

class TetViewerState():
    """
    storage for the current viewing state
    """
    def __init__(self):
        """
        initalize the state
        """
        self._surface_ctr = None
        self._current_tet_ctr = None
        self._euler_x = 0.0
        self._euler_y = 0.0
        self._clear_colour = None
        self._ctr_is_tet = True
        self._shift = (0.0, 0.0, 0.0)

    def reset(self):
        """
        reset the shift and rotation
        """
        self._euler_x = 0.0
        self._euler_y = 0.0
        self._shift = (0.0, 0.0, 0.0)

    def get_clear_color(self):
        """
        get the current clear color
        Returns:
            tuple(float): (r, g, b) 0.0 <= number <=1.0
        """
        return self._clear_colour

    def set_clear_color(self, red, green, blue):
        """
        set the clear color
        Args:
            red: float 0.0<=red<=1.0
            green: float 0.0<=green<=1.0
            blue: float 0.0<=blue<=1.0
        """
        self._clear_colour = (red, green, blue)

    def get_current_ctr(self):
        """
        get the currently centred point
        Returns:
            tuple(float): (x, y, z)
        """
        if self._ctr_is_tet:
            return self._current_tet_ctr

        return self._surface_ctr

    def current_ctr_is_tet(self):
        """
        is the current ctr a tet
        Return:
            (bool): True if current ctr is a tet else False
        """
        return self._ctr_is_tet

    def center_on_tet(self):
        """
        set the current ctr to be a tet
        """
        self._ctr_is_tet = True

    def center_on_surface(self):
        """
        set the current ctr to be a surface
        """
        self._ctr_is_tet = False

    def set_tet_ctr(self, x, y, z):
        """
        set the centre of the current tet
        Args:
            x: float
            y: float
            z: float
        """
        self._current_tet_ctr = (x, y, z)

    def set_surface_ctr(self, x, y, z):
        """
        set the centre of the current tet
        Args:
            x: float
            y: float
            z: float
        """
        self._surface_ctr = (x, y, z)

    def get_euler_x(self):
        """
        get the current rotation about x axis
        Returns
            float: rotation about x in degrees
        """
        return self._euler_x

    def get_euler_y(self):
        """
        get the current rotation about y axis
        Returns
            float: rotation about y in degrees
        """
        return self._euler_y

    def set_euler_x(self, angle):
        """
        set rotation about x axis
        Args:
            angle: float in degrees
        """
        self._euler_x = angle

    def set_euler_y(self, angle):
        """
        set rotation about y axis
        Args:
            angle: float in degrees
        """
        self._euler_y = angle

    def get_shift(self):
        """
        get the current shift
        Returns
            tuple float: (x, y, z)
        """
        return self._shift

    def set_shift_xy(self, x, y):
        """
        set xy componants of shift
        Args:
            x: float
            y: float
        """
        self._shift = (x, y, self._shift[2])

    def set_shift_z(self, z):
        """
        set z componants of shift
        Args:
            z: float
        """
        self._shift = (self._shift[0], self._shift[1], z)

    def set_shift(self, x, y, z):
        """
        set the shift
        Args:
            x: float
            y: float
            z: float
        """
        self._shift = (x, y, z)
