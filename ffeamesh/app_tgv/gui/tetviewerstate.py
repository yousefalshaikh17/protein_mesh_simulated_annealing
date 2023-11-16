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
# set up linting conditions
# pylint: disable = import-error
import json

import ffeamesh.tetprops as tp

## z for centring
CTR_Z = 1350.0

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
        self._current_tet_nodes = None
        self._display_current_tet = True
        self._current_ctr = None

        self._euler_x = 0.0
        self._euler_y = 0.0
        self._clear_colour = None
        self._shift = (0.0, 0.0, 0.0)
        self._look_z = -1500.0
        self._edges_width = 1
        self._perspective = True

        ## highlight peak nodes
        self._show_peak_nodes = False

    def reset(self):
        """
        reset the shift and rotation
        """
        self._euler_x = 0.0
        self._euler_y = 0.0
        self._shift = (0.0, 0.0, 0.0)

    def set_edges_width(self, value):
        """
        setter for the screen width of the lines in the edges
        Args:
            value (int): line width in pixels
        """
        self._edges_width = value

    def get_edges_width(self):
        """
        getter for the screen width of the lines in the edges
        Returns:
            int: line width in pixels
        """
        return self._edges_width

    def set_perspective(self, flag):
        """
        set the perspective flag
        Args:
            flag: bool
        """
        if not isinstance(flag, bool):
            raise ValueError("attempt to set TetViewState._perspectiv to non-bool")

        self._perspective = flag

    def perspective_view(self):
        """
        get the perspective flag if True use perspective
        Args:
            flag: bool
        """
        return self._perspective

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
        return self._current_ctr

    def center_on_tet(self):
        """
        set the current ctr to be a tet
        """
        self._current_ctr = self._current_tet_ctr
        self._shift = (self._current_tet_ctr[0],
                       self._current_tet_ctr[1],
                       CTR_Z)

    def center_on_surface(self):
        """
        set the current ctr to be a surface
        """
        self._current_ctr = self._surface_ctr
        self._shift = (self._surface_ctr[0],
                       self._surface_ctr[1],
                       CTR_Z)

    def set_current_tet(self, tet):
        """
        set the current tet
        """
        self._current_tet_nodes = tet
        tmp = tp.find_centre(tet)
        self._current_tet_ctr = (tmp[0], tmp[1], tmp[2])

    def clear_current_tet(self):
        """
        remove the current tet
        """
        self._current_tet_nodes = None
        self._current_tet_ctr = None

    def get_current_tet(self):
        """
        get the current tet
        """
        return self._current_tet_nodes

    def set_surface_ctr(self, x_coord, y_coord, z_coord):
        """
        set the centre of the current tet
        Args:
            x_coord: float
            y_coord: float
            z_coord: float
        """
        self._surface_ctr = (x_coord, y_coord, z_coord)
        self.center_on_surface()

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

    def set_shift_xy(self, x_coord, y_coord):
        """
        set xy componants of shift
        Args:
            x_coord: float
            y_coord: float
        """
        self._shift = (x_coord, y_coord, self._shift[2])

    def set_shift_z(self, z_coord):
        """
        set z componants of shift
        Args:
            z_coord: float
        """
        self._shift = (self._shift[0], self._shift[1], z_coord)

    def set_shift(self, x_coord, y_coord, z_coord):
        """
        set the shift
        Args:
            x_coord: float
            y_coord: float
            z_coord: float
        """
        self._shift = (x_coord, y_coord, z_coord)

    def get_look_z(self):
        """
        get the z componant of the look-from
        Returns
            float
        """
        return self._look_z

    def set_look_z(self, z_coord):
        """
        set the z componant of the look-from
        Args:
            z_coord: float
        """
        self._look_z = z_coord

    def set_display_current_tet(self, flag):
        """
        set the tet to be displayed
        """
        self._display_current_tet = flag

    def display_current_tet(self):
        """
        display current tet, if there is one
        """
        if self._display_current_tet and self._current_tet_nodes is not None:
            return True

        return False

    def toggle_show_peak_nodes(self):
        """
        flip the show peak nodes flag
        """
        self._show_peak_nodes = not self._show_peak_nodes

    def get_show_peak_nodes(self):
        """
        get the show peak nodes flag
        """
        return self._show_peak_nodes

    def save_setup(self, file_path):
        """
        save the view
        Args:
            file_path (pathlib.Path)
        """
        state = json.dumps(self, cls=EncodeViewState)
        with file_path.open('w', encoding='utf-8') as file:
            print(state, file=file)

    def load_setup(self, file_path):
        """
        load an object from a JSON string
        Args:
            obj (TetViewerState): object to be loded
            file_path (pathlib.Path): file holding json data
        """
        data = None
        with file_path.open('r', encoding='utf-8') as file:
            data = json.load(file)

        shift = data["shift"]
        self.set_shift(shift[0], shift[1], shift[2])
        self.set_euler_x(data["euler_x"])
        self.set_euler_y(data["euler_y"])

class EncodeViewState(json.JSONEncoder):
    """
    json encoding object
    """
    def default(self, o):
        """
        encoding function converts object to dict
        Args:
            o (TetViewerState): object to be encoded
        Returns:
            dict: holding json data
        Raises:
            TypeError if obj is not TetViewerState
        """
        if not isinstance(o, TetViewerState):
            raise TypeError(f"object is not viewstate.ViewState: was {type(o)}")

        return {'shift': o.get_shift(),
                'euler_x': o.get_euler_x(),
                'euler_y': o.get_euler_y()}
