"""
subclass of QOpenGLWidget (drawing area) provides all drawing methods

----------------------------------------------

Licensed under the GNU General Public License, Version 3.0 (the "License"); you
may not use this file except in compliance with the License. You may obtain a
copy of the License at <https://www.gnu.org/licenses/gpl-3.0.html>.

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
# pylint: disable = invalid-name

from enum import Enum
import numpy as np

import PyQt5.QtWidgets as qw
import PyQt5.QtCore as qc
import PyQt5.QtGui as qg
import PyQt5.Qt as qt

import OpenGL.GL as gl
import OpenGL.GLU as glu

from tetmeshtools.app_tgv.gui.sphere import Sphere
from tetmeshtools.app_tgv.gui.tetviewerstate import TetViewerState

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
    subclass of QOpenGLWidget (drawing area / graphics context) provides all drawing methods
    """

    ## notify that rotation user input needs to be reset
    reset_input = qc.pyqtSignal()

    def __init__(self, parent=None):
        """initalize the window"""
        self.parent = parent
        super().__init__(parent)

        ## the state
        self._state = TetViewerState()

        ## the sphere
        self._sphere = Sphere()

        ## storage for show surfaces
        self._model = None

        ## storage for the mouse button states (used in dragging & zooming)
        self._mouse_state    = MouseStates.NONE

        ## storae for the current mouse position (used in dragging & zooming)
        self._mouse_position = None

        ## the minimum z value allowed in zooming
        self._min_z = -2000.0

        ## the maximum z value allowd in zooming
        self._max_z = 2000.0

        ## z coord of the near clipping plane
        self._near_clip = 1.0

        ## z coord of the far clipping plane
        self._far_clip = 3000.0

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
            glu.gluPerspective(self._state.get_field_of_view(),
                               aspect,
                               self._near_clip,
                               self._far_clip)
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

        # OpenGL.constant.IntConstant not regarded by pylint as an int hence warning
        gl.glClear(gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT)

        if self._model is None:
            return

        gl.glPushMatrix()    # push the current matrix to the current stack

        glu.gluLookAt(0.0, 0.0, self._state.get_look_from_z(),
                      0.0, 0.0, 0.0,
                      0.0, 1.0, 0.0)

        # apply the user's current shift (up, down, zoom)
        shift = self._state.get_shift()
        gl.glTranslate(-shift[0], -shift[1], -shift[2])

        # apply the current rotations (origin in model)
        gl.glRotate(self._state.get_euler_x(), 1.0, 0.0, 0.0)
        gl.glRotate(self._state.get_euler_y(), 0.0, 1.0, 0.0)

        # move mesh ctr to origin
        ctr = self._state.get_current_ctr()
        gl.glTranslate(-ctr[0], -ctr[1], -ctr[2])

        if self._state.display_current_tet():
            self.draw_selected_tet()
        if self._state.get_show_lattice():
            self.draw_triangle_outline()
        if self._state.get_show_faces():
            self.draw_triangles()

        gl.glPopMatrix()    # restore the previous modelview matrix

    def draw_triangles(self, scale=(1.0, 1.0, 1.0)):
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

    def draw_triangle_outline(self, scale=(1.0, 1.0, 1.0)):
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

    def draw_selected_tet(self, scale=(1.0, 1.0, 1.0)):
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

    @qc.pyqtSlot(int)
    def show_faces(self, check_state):
        """
        show the faces
        Args:
            check_state (Qt.CheckState)
        """
        if check_state == qt.Qt.CheckState.Checked:
            self._state.set_show_faces(True)
        else:
            self._state.set_show_faces(False)

        self.update()

    @qc.pyqtSlot(int)
    def show_surface_lattice(self, check_state):
        """
        show the faces edges
        Args:
            check_state (Qt.CheckState)
        """
        if check_state == qt.Qt.CheckState.Checked:
            self._state.set_show_lattice(True)
        else:
            self._state.set_show_lattice(False)

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
        new_shift_z = shift[2] - float(del_y)
        effective_z = self._state.get_look_from_z() + new_shift_z

        self._mouse_position = mouse_position

        if self._max_z > effective_z > self._min_z:
            self._state.set_shift_z(new_shift_z)

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
        set the view matrix (persepctive, orthogonal)
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

        ctr, radius = self._model.get_bounding_sphere()
        radius = radius*1.1

        view_distance = radius/np.tan(np.radians(self._state.get_field_of_view()/2.0))

        self._max_z = ctr[2] - radius
        self._min_z = ctr[2] - self._far_clip

        self._state.set_look_from_z(ctr[2]-view_distance)
        self._state.set_surface_ctr(ctr[0], ctr[1], ctr[2])
        self._state.centre_on_surface()
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

    @qc.pyqtSlot(float, float)
    def set_min_max_z(self, min_z, max_z):
        """
        set the minimum and maximum limits on the z axis
        Args:
            min_z (float): minimum allowed z value
            max_z (float): maximum allowed z value
        """
        self._max_z = min_z
        self._max_z = max_z

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
            self._state.centre_on_surface()
            self.update()

    @qc.pyqtSlot()
    def centre_tet(self):
        """
        callback for click of ctr tet button
        """
        if self._model is not None and self._state.get_current_tet() is not None:
            self._state.centre_on_tet()
            self.update()

    @qc.pyqtSlot(int)
    def show_current_tet(self, check_state):
        """
        tobggle display of current tet
        Args:
            check_state (Qt.CheckState)
        """
        if check_state == qt.Qt.CheckState.Unchecked:
            self._state.set_display_current_tet(False)
        else:
            self._state.set_display_current_tet(True)
        self.update()
