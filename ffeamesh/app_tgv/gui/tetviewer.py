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
# set up linting conditions
# pylint: disable = import-error
# pylint: disable = c-extension-no-member
import numpy as np
from enum import Enum

import PyQt5.QtWidgets as qw
import PyQt5.QtCore as qc
import PyQt5.QtGui as qg

import OpenGL.GL as gl
import OpenGL.GLU as glu

from ffeamesh.app_tgv.gui.sphere import Sphere
from ffeamesh.app_tgv.gui.tetviewerstate import TetViewerState

#TODO get rid of 3 state variable faces lattice tet and use call to parent in paint
#     replace Hide tet button with check box in parent

class MouseStates(Enum):
    """
    possible states of the mouse input
    """
    ## not important
    NONE = 0

    ## motion
    MOTION = 1

    ## zooming
    ZOOM = 2

class TetViewer(qw.QOpenGLWidget):
    """
    the view on graphics context
    """

    ## notify that rotation user input needs to be reset
    reset_input = qc.pyqtSignal()

    def __init__(self, parent=None):
        """initalize the window"""
        self.parent = parent
        qw.QOpenGLWidget.__init__(self, parent)

        ## the state
        self._state = TetViewerState()

        ## the sphere
        self._sphere = Sphere()

        ## storage for show surfaces
        self._model = None

        ## show/hide surface faces
        self._show_faces = False

        ## show/hide surface lattice
        self._show_lattice = False

        self._mouse_state    = MouseStates.NONE
        self._mouse_position = None

        self.set_background()

    def initializeGL(self):
        """initalize the graphics context"""
        gl.glEnable(gl.GL_DEPTH_TEST)               # enable OGL depth testing

    def resizeGL(self, width, height):
        """
        opengl resize callback
        Arts
            width (int): width in pixels
            height (int): height in pixels
        """
        gl.glViewport(0, 0, width, height)
        self.set_projection(width, height)

    def set_projection(self, width, height):
        """
        set the projection matrix
        Arts
            width (int): width in pixels
            height (int): height in pixels
        """
        gl.glMatrixMode(gl.GL_PROJECTION)
        gl.glLoadIdentity()
        aspect = width / float(height)

        if self._state.perspective_view():
            glu.gluPerspective(45.0, aspect, 1.0, 3000.0)
        else:
            gl.glOrtho(-500.0, 500.0, -500.0, 500.0, 15.0, 3000.0)

        gl.glMatrixMode(gl.GL_MODELVIEW)
        gl.glLoadIdentity()

    def paintGL(self):
        """
        render the contents callback
        """
        c_color = self._state.get_clear_color()
        gl.glClearColor(c_color[0], c_color[1], c_color[2], 1.0)
        gl.glClear(gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT)

        if self._model is None:
            return

        gl.glPushMatrix()    # push the current matrix to the current stack

        glu.gluLookAt(0.0, 0.0, self._state.get_look_z(),
                      0.0, 0.0, 0.0,
                      0.0, 1.0, 0.0)

        # location of tet as input by user
        shift = self._state.get_shift()
        gl.glTranslate(-shift[0], -shift[1], -shift[2])

        ctr = self._state.get_current_ctr()
        gl.glTranslate(ctr[0], ctr[1], ctr[2])
        gl.glRotate(self._state.get_euler_x(), 1.0, 0.0, 0.0)
        gl.glRotate(self._state.get_euler_y(), 0.0, 1.0, 0.0)
        gl.glTranslate(-ctr[0], -ctr[1], -ctr[2])

        scale = (1.0, 1.0, 1.0)
        if self._state.display_current_tet():
            self.draw_selected_tet(scale)
        if self._show_lattice:
            self.draw_triangle_outline(scale)
        if self._show_faces:
            self.draw_triangles(scale)

        gl.glPopMatrix()    # restore the previous modelview matrix

    def draw_triangles(self, scale):
        """
        render the triangles
        Args:
            scale (tuple float)
        """
        gl.glEnable(gl.GL_BLEND)
        gl.glBlendFunc(gl.GL_ONE, gl.GL_SRC_ALPHA)
        gl.glBlendEquation(gl.GL_FUNC_ADD)

        gl.glPushMatrix()
        gl.glScale(scale[0], scale[1], scale[2])
        gl.glPushAttrib(gl.GL_COLOR_BUFFER_BIT)

        mat_specular = [1.0, 1.0, 1.0, 1.0]
        mat_shininess = [50.0]
        mat_amb_diff = [0.1, 0.5, 0.8, 0.5]

        light_ambient = [0.0, 0.0, 0.0, 1.0]
        light_diffuse = [1.0, 1.0, 1.0, 1.0]
        light_specular = [1.0, 1.0, 1.0, 1.0]

        light_position = [1000.0, 100.0, 100.0, 0.0]
        gl.glShadeModel (gl.GL_SMOOTH)
        gl.glMaterialfv(gl.GL_FRONT_AND_BACK, gl.GL_AMBIENT_AND_DIFFUSE, mat_amb_diff)
        gl.glMaterialfv(gl.GL_FRONT, gl.GL_SPECULAR, mat_specular)
        gl.glMaterialfv(gl.GL_FRONT, gl.GL_SHININESS, mat_shininess)

        gl.glLightfv(gl.GL_LIGHT0, gl.GL_AMBIENT, light_ambient)
        gl.glLightfv(gl.GL_LIGHT0, gl.GL_DIFFUSE, light_diffuse)
        gl.glLightfv(gl.GL_LIGHT0, gl.GL_SPECULAR, light_specular)
        gl.glLightfv(gl.GL_LIGHT0, gl.GL_POSITION, light_position)

        gl.glEnable(gl.GL_LIGHTING)
        gl.glEnable(gl.GL_LIGHT0)
        gl.glBegin(gl.GL_TRIANGLES)

        surface = self._model.get_surface()
        faces = surface.get_faces()
        for index in faces:
            nodes = surface.get_triangle_nodes(index)

            edge0 = nodes[0].to_edge_array(nodes[1])
            edge1 = nodes[1].to_edge_array(nodes[2])
            normal = np.cross(edge0, edge1)
            norm = np.linalg.norm(normal)
            normal = normal/norm

            gl.glNormal(normal[0], normal[1], normal[2])
            gl.glVertex(nodes[0].x, nodes[0].y, nodes[0].z)
            gl.glVertex(nodes[1].x, nodes[1].y, nodes[1].z)
            gl.glVertex(nodes[2].x, nodes[2].y, nodes[2].z)

        gl.glEnd()

        gl.glDisable(gl.GL_LIGHT0)
        gl.glDisable(gl.GL_LIGHTING)

        gl.glPopAttrib()
        gl.glPopMatrix()
        gl.glDisable(gl.GL_BLEND)

    def draw_triangle_outline(self, scale):
        """
        render the triangles outlines
        Args:
            scale (tuple float)
        """
        gl.glPushMatrix()
        gl.glScale(scale[0], scale[1], scale[2])
        gl.glPushAttrib(gl.GL_COLOR_BUFFER_BIT)

        gl.glColor(0.1, 0.0, 0.7, 1.0)
        gl.glLineWidth(self._state.get_edges_width())
        gl.glBegin(gl.GL_LINES)
        surface = self._model.get_surface()
        faces = surface.get_faces()
        for index in faces:
            nodes = surface.get_triangle_nodes(index)
            gl.glVertex3f(nodes[0].x, nodes[0].y, nodes[0].z)
            gl.glVertex3f(nodes[1].x, nodes[1].y, nodes[1].z)
            gl.glVertex(nodes[1].x, nodes[1].y, nodes[1].z)
            gl.glVertex(nodes[2].x, nodes[2].y, nodes[2].z)
            gl.glVertex(nodes[2].x, nodes[2].y, nodes[2].z)
            gl.glVertex(nodes[0].x, nodes[0].y, nodes[0].z)
        gl.glEnd()

        gl.glPopAttrib()
        gl.glPopMatrix()

    def draw_selected_tet(self, scale):
        """
        draw the user selecte tet
        """
        current_tet = self._state.get_current_tet()

        self._sphere.draw(current_tet[0], scale)
        self._sphere.draw(current_tet[1], scale)
        self._sphere.draw(current_tet[2], scale)
        self._sphere.draw(current_tet[3], scale)

        gl.glLineWidth(10.0)
        gl.glBegin(gl.GL_LINES)
        gl.glColor3f(0.0, 1.0, 0.0)
        gl.glVertex3f(current_tet[0].x, current_tet[0].y, current_tet[0].z)
        gl.glVertex3f(current_tet[1].x, current_tet[1].y, current_tet[1].z)

        gl.glVertex3f(current_tet[0].x, current_tet[0].y, current_tet[0].z)
        gl.glVertex3f(current_tet[2].x, current_tet[2].y, current_tet[2].z)

        gl.glVertex3f(current_tet[0].x, current_tet[0].y, current_tet[0].z)
        gl.glVertex3f(current_tet[3].x, current_tet[3].y, current_tet[3].z)

        gl.glVertex3f(current_tet[1].x, current_tet[1].y, current_tet[1].z)
        gl.glVertex3f(current_tet[2].x, current_tet[2].y, current_tet[2].z)

        gl.glVertex3f(current_tet[1].x, current_tet[1].y, current_tet[1].z)
        gl.glVertex3f(current_tet[3].x, current_tet[3].y, current_tet[3].z)

        gl.glVertex3f(current_tet[2].x, current_tet[2].y, current_tet[2].z)
        gl.glVertex3f(current_tet[3].x, current_tet[3].y, current_tet[3].z)

        gl.glEnd()

    def display_tet(self, tet_index):
        """
        display a tet
        Args:
            tet_indes (int): the tet's index
        """
        self._state.set_current_tet(self._model.get_tet_nodes(tet_index))
        self.reset_view()
        self.update()

    def reset_all(self):
        """
        reset back to empty
        """
        self._state = TetViewerState()
        self._model = None
        self.set_background()

        if self.isValid():
            self.makeCurrent()
            self.set_projection(self.width(), self.height())
            self.update()

    def show_faces(self, flag):
        """
        show the faces
        Args:
            flag (bool)
        """
        self._show_faces = flag
        self.update()

    def show_surface_lattice(self, flag):
        """
        show the faces edges
        Args:
            flag (bool)
        """
        self._show_lattice = flag
        self.update()

    def reset_view(self):
        """
        reset view parameters
        """
        self._state.reset()

    def save_setup(self, file_path):
        """
        save the view setup
        Args:
            file_path (pathlib.Path)
        """
        self._state.save_setup(file_path)

    def load_setup(self, file_path):
        """
        load a view setup
        Args:
            file_path (pathlib.Path)
        """
        self._state.load_setup(file_path)
        self.update()

    def shift(self, mouse_position):
        """
        move camera 'up-down'
        """
        old_shift = self._state.get_shift()
        del_x = self._mouse_position.x() - mouse_position.x()
        del_y = mouse_position.y() - self._mouse_position.y()
        new_x = old_shift[0] - float(del_x)/5.0
        new_y = old_shift[1] + float(del_y)/5.0

        self._state.set_shift_xy(new_x, new_y)

        self._mouse_position = mouse_position

        self.update()

    def zoom(self, mouse_position):
        """
        move camera 'in-out'
        """
        if not self._state.perspective_view():
            return

        del_y = mouse_position.y() - self._mouse_position.y()

        shift = self._state.get_shift()
        new_z = shift[2] - float(del_y)

        self._mouse_position = mouse_position
        if new_z > -2000.0:
            self._state.set_shift_z(new_z)

        self.update()

    def set_background(self, text="White"):
        """
        set the background colour
        """
        if text == "Black":
            self._state.set_clear_color(0.0, 0.0, 0.0)
        elif text == "White":
            self._state.set_clear_color(1.0, 1.0, 1.0)
        elif text == "Gray":
            self._state.set_clear_color(0.3, 0.4, 0.5)

    def change_background(self, text):
        """
        change the clear color
        Args:
            text (str): one of {"Black", "White", "Gray"}
        """
        self.set_background(text)
        self.update()

    def set_view(self, text):
        """
        set the background colour
        """
        if text == "Perspective":
            self._state.set_perspective(True)
        else:
            self._state.set_perspective(False)

        self.makeCurrent()
        self.set_projection(self.width(), self.height())

        self.update()

    def set_model(self, model):
        """
        set the current model
        Args:
            model (TetModel)
        """
        self._model = model
        ctr = self._model.get_surface().get_surface_ctr()
        self._state.set_surface_ctr(ctr[0], ctr[1], ctr[2])
        self._state.center_on_surface()
        self._state.clear_current_tet()

    @qc.pyqtSlot(int)
    def set_thickness(self, val):
        """
        set thickness of surface triangle edges
        """
        self._state.set_edges_width(val)
        self.update()

    @qc.pyqtSlot(int)
    def set_rot_x(self, val):
        """
        set rotation about x axis
        """
        self._state.set_euler_x(float(val))
        self.update()

    @qc.pyqtSlot(int)
    def set_rot_y(self, val):
        """
        set rotation about y axis
        """
        self._state.set_euler_y(float(val))
        self.update()

    @qc.pyqtSlot(qg.QMouseEvent)
    def mouseMoveEvent(self, event):
        """
        callback for the motion of a mouse with a button held
        """
        if self._mouse_state == MouseStates.ZOOM:
            self.zoom(event.pos())
        elif self._mouse_state == MouseStates.MOTION:
            self.shift(event.pos())

        return super().mouseMoveEvent(event)

    @qc.pyqtSlot(qg.QMouseEvent)
    def mousePressEvent(self, event):
        """
        callback for the press of a mouse button
        """
        if event.buttons() in (qc.Qt.LeftButton, qc.Qt.RightButton):
            self._mouse_position = event.pos()
            if event.buttons() == qc.Qt.LeftButton:
                self._mouse_state = MouseStates.MOTION
            else:
                self._mouse_state = MouseStates.ZOOM

        return super().mousePressEvent(event)

    @qc.pyqtSlot(qg.QMouseEvent)
    def mouseReleaseEvent(self, event):
        """
        callback for the release of a mouse button
        """
        self._mouse_state = MouseStates.NONE
        self._mouse_position = None

        return super().mouseReleaseEvent(event)

    @qc.pyqtSlot()
    def centre_mesh(self):
        """
        callback for click of ctr mesh button
        """
        if self._model is not None:
            self._state.center_on_surface()
            self.update()

    @qc.pyqtSlot()
    def centre_tet(self):
        """
        callback for click of ctr tet button
        """
        if self._model is not None and self._state.get_current_tet() is not None:
            self._state.center_on_tet()
            self.update()

    @qc.pyqtSlot(bool)
    def show_current_tet(self, flag):
        """
        stop displaying current tet
        """
        self._state.set_display_current_tet(flag)
        self.update()
