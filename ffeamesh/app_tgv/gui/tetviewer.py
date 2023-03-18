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
from enum import Enum

import PyQt5.QtWidgets as qw
import PyQt5.QtCore as qc
import PyQt5.QtGui as qg

import OpenGL.GL as gl
import OpenGL.GLU as glu

from ffeamesh.app_tgv.gui.sphere import Sphere
from ffeamesh.app_tgv.gui.threestore import ThreeStore

import ffeamesh.tetprops as tp

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
    reset_rot_input = qc.pyqtSignal()

    def __init__(self, parent=None):
        """initalize the window"""
        self.parent = parent
        qw.QOpenGLWidget.__init__(self, parent)

        ## euler angles
        self._rotations =ThreeStore(0.0, 0.0, 0.0)

        ## look from point
        self._look_z = -1500.0

        ## perspective view
        self._perspective = True

        self._shift = ThreeStore(0.0, 0.0, 0.0)

        ## the sphere
        self._sphere = Sphere()

        ## storage for show surfaces
        self._surface = None

        ## storage for show lattice
        self._lattice = None

        ## flag to indicat a reset event
        self._reset = False

        ## the four vertices of the tet being viewed
        self._current_tet = None

        ## centre of the tet being viewed
        self._current_tet_ctr = None

        self._mouse_state    = MouseStates.NONE
        self._mouse_position = None

        self._clear_colour = None
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

        if self._perspective:
            glu.gluPerspective(45.0, aspect, 1.0, 3000.0)
        else:
            gl.glOrtho(-500.0, 500.0, -500.0, 500.0, 15.0, 3000.0)

        gl.glMatrixMode(gl.GL_MODELVIEW)
        gl.glLoadIdentity()

    def paintGL(self):
        """
        render the contents callback
        """
        gl.glClearColor(self._clear_colour.x, self._clear_colour.y, self._clear_colour.z, 1)
        gl.glClear(gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT)

        gl.glPushMatrix()    # push the current matrix to the current stack

        glu.gluLookAt(0.0, 0.0, self._look_z,
                      0.0, 0.0, 0.0,
                      0.0, 1.0, 0.0)

        # location of tet as input by user
        gl.glTranslate(-self._shift.x, -self._shift.y, self._shift.z)

        if self._current_tet is not None:
            gl.glTranslate(self._current_tet_ctr.x,
                           self._current_tet_ctr.y,
                           self._current_tet_ctr.z)
            gl.glRotate(self._rotations.x, 1.0, 0.0, 0.0)
            gl.glRotate(self._rotations.y, 0.0, 1.0, 0.0)
            gl.glTranslate(-self._current_tet_ctr.x,
                           -self._current_tet_ctr.y,
                           -self._current_tet_ctr.z)

        scale = ThreeStore(1.0, 1.0, 1.0)
        if self._current_tet is not None:
            self.draw_selected_tet(scale)
        if self._lattice is not None:
            self.draw_triangle_outline(scale)
        if self._surface is not None:
            self.draw_triangles(scale)

        gl.glPopMatrix()    # restore the previous modelview matrix

    def draw_triangles(self, scale):
        """
        render the triangles
        """
        import numpy as np
        nodes = self._surface["nodes"]
        gl.glEnable(gl.GL_BLEND)
        gl.glBlendFunc(gl.GL_ONE, gl.GL_SRC_ALPHA)
        gl.glBlendEquation(gl.GL_FUNC_ADD)

        gl.glPushMatrix()
        gl.glScale(scale.x, scale.y, scale.z)
        gl.glPushAttrib(gl.GL_COLOR_BUFFER_BIT)

        #gl.glColor4f(1.0, 0.0, 0.7, 0.5)
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
        for index, verts_indices in self._surface.items():
            if index != "nodes":
                node0 = nodes[verts_indices[0]]
                node1 = nodes[verts_indices[1]]
                node2 = nodes[verts_indices[2]]

                edge0 = node0.to_edge_array(node1)
                edge1 = node1.to_edge_array(node2)
                normal = np.cross(edge0, edge1)
                norm = np.linalg.norm(normal)
                normal = normal/norm

                gl.glNormal(normal[0], normal[1], normal[2])
                gl.glVertex(node0.x, node0.y, node0.z)
                gl.glVertex(node1.x, node1.y, node1.z)
                gl.glVertex(node2.x, node2.y, node2.z)


        gl.glEnd()

        gl.glDisable(gl.GL_LIGHT0)
        gl.glDisable(gl.GL_LIGHTING)

        gl.glPopAttrib()
        gl.glPopMatrix()
        gl.glDisable(gl.GL_BLEND)

    def draw_triangle_outline(self, scale):
        """
        render the triangles outlines
        """
        nodes = self._lattice["nodes"]

        gl.glPushMatrix()
        gl.glScale(scale.x, scale.y, scale.z)
        gl.glPushAttrib(gl.GL_COLOR_BUFFER_BIT)

        gl.glColor(1.0, 0.0, 0.0, 1.0)
        gl.glLineWidth(5)
        gl.glBegin(gl.GL_LINES)
        for index, verts_indices in self._lattice.items():
            if index != "nodes":
                node0 = nodes[verts_indices[0]]
                node1 = nodes[verts_indices[1]]
                node2 = nodes[verts_indices[2]]
                gl.glVertex3f(node0.x, node0.y, node0.z)
                gl.glVertex3f(node1.x, node1.y, node1.z)
                gl.glVertex(node1.x, node1.y, node1.z)
                gl.glVertex(node2.x, node2.y, node2.z)
                gl.glVertex(node2.x, node2.y, node2.z)
                gl.glVertex(node0.x, node0.y, node0.z)
        gl.glEnd()

        gl.glPopAttrib()
        gl.glPopMatrix()

    def draw_selected_tet(self, scale):
        """
        draw the user selecte tet
        """
        one = self._current_tet[0]
        two = self._current_tet[1]
        three = self._current_tet[2]
        four = self._current_tet[3]

        self._sphere.draw(one, scale)
        self._sphere.draw(two, scale)
        self._sphere.draw(three, scale)
        self._sphere.draw(four, scale)

        gl.glLineWidth(10.0)
        gl.glBegin(gl.GL_LINES)
        gl.glColor3f(0.0, 1.0, 0.0)
        gl.glVertex3f(one.x, one.y, one.z)
        gl.glVertex3f(two.x, two.y, two.z)

        gl.glVertex3f(one.x, one.y, one.z)
        gl.glVertex3f(three.x, three.y, three.z)

        gl.glVertex3f(one.x, one.y, one.z)
        gl.glVertex3f(four.x, four.y, four.z)

        gl.glVertex3f(two.x, two.y, two.z)
        gl.glVertex3f(three.x, three.y, three.z)

        gl.glVertex3f(two.x, two.y, two.z)
        gl.glVertex3f(four.x, four.y, four.z)

        gl.glVertex3f(three.x, three.y, three.z)
        gl.glVertex3f(four.x, four.y, four.z)

        gl.glEnd()

    def display(self, tet, nodes):
        """
        display a tet
        Args:
            tet (tetgen_read.Tetrahedron4)
            nodes ({NodePoint}): the points
        """
        self._current_tet = []
        self._current_tet.append(nodes[tet.vert0])
        self._current_tet.append(nodes[tet.vert1])
        self._current_tet.append(nodes[tet.vert2])
        self._current_tet.append(nodes[tet.vert3])

        ctr = tp.find_centre(self._current_tet)
        self._current_tet_ctr = ThreeStore(ctr[0], ctr[1], ctr[2])

        self.reset_view()
        self.update()

    def reset_all(self):
        """
        reset back to empty
        """
        self._current_tet = None
        self._current_tet_ctr = None
        self._surface = None
        self._lattice = None
        self._clear_colour = None
        self._perspective = True
        self.set_background()

        if self.isValid():
            self.makeCurrent()
            self.set_projection(self.width(), self.height())
            self.update()

    def show_faces(self, nodes, faces):
        """
        show the faces
        Args:
            nodes
            faces
        """
        self._surface = {}
        self._surface["nodes"] = nodes
        for face in faces:
            self._surface[face.index] = [face.vert0, face.vert1, face.vert2]

        self.update()

    def show_surface_lattice(self, nodes, faces):
        """
        show the faces
        Args:
            nodes
            faces
        """
        self._lattice = {}
        self._lattice["nodes"] = nodes
        for face in faces:
            self._lattice[face.index] = [face.vert0, face.vert1, face.vert2]

        self.update()

    def hide_faces(self):
        """
        stop showing faces
        """
        self._surface = None
        self.update()

    def hide_surface_lattice(self):
        """
        stop showing faces
        """
        self._lattice = None
        self.update()

    def reset_view(self):
        """
        reset view parameters
        """
        self._rotations.x = 0.0
        self._rotations.y = 0.0
        self._rotations.z = 0.0

        self._reset = True
        self.reset_rot_input.emit()

    def shift(self, mouse_position):
        """
        move camera 'up-down'
        """
        del_x = self._mouse_position.x() - mouse_position.x()
        del_y = mouse_position.y() - self._mouse_position.y()

        self._shift.x -= float(del_x)/5.0
        self._shift.y += float(del_y)/5.0

        self._mouse_position = mouse_position

        self.update()

    def zoom(self, mouse_position):
        """
        move camera 'in-out'
        """
        if not self._perspective:
            return

        del_y = mouse_position.y() - self._mouse_position.y()

        new_z = self._shift.z - float(del_y)

        self._mouse_position = mouse_position
        if new_z > -2000.0: #-1490.0:
            self._shift.z = new_z

        self.update()

    def set_background(self, text="White"):
        """
        set the background colour
        """
        if text == "Black":
            self._clear_colour = ThreeStore(0.0, 0.0, 0.0)
        elif text == "White":
            self._clear_colour = ThreeStore(1.0, 1.0, 1.0)
        elif text == "Gray":
            self._clear_colour = ThreeStore(0.3, 0.4, 0.5)

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
            self._perspective = True
        else:
            self._perspective = False

        self.makeCurrent()
        self.set_projection(self.width(), self.height())

        self.update()

    @qc.pyqtSlot(int)
    def set_rot_x(self, val):
        """
        set rotation about x axis
        """
        self._rotations.x = int(val)
        self.update()

    @qc.pyqtSlot(int)
    def set_rot_y(self, val):
        """
        set rotation about y axis
        """
        self._rotations.y = int(val)
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
    def centre_tet(self):
        """
        callback for click of ctr button
        """
        if self._current_tet_ctr is not None:
            self._shift.x = self._current_tet_ctr.x
            self._shift.y = self._current_tet_ctr.y
            self._shift.z = -1350.0

        self.update()

    @qc.pyqtSlot()
    def remove_current_tet(self):
        """
        stop displaying current tet
        """
        if self._current_tet_ctr is not None:
            self._current_tet_ctr = None
            self._current_tet     = None
            self.reset_view()
            self.update()
