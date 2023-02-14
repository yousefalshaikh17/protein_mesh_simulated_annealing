"""
Created on 03 Jan 2023

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
# pylint: disable = c-extension-no-member
import OpenGL.GL as gl
import OpenGL.GLU as glu

class Sphere():
    """a graphics sphere"""
    def __init__(self):
        """initialize the object"""
        self._quadric = glu.gluNewQuadric()

    def draw(self, translate, scale, color=None):
        """
        render the sphere
        """
        if color is None:
            color = [0.0, 0.1, 0.9, 0.1]
        gl.glPushMatrix()
        gl.glTranslate(translate.x, translate.y, translate.z)
        gl.glScale(scale.x, scale.y, scale.z)
        gl.glPushAttrib(gl.GL_COLOR_BUFFER_BIT)
        gl.glEnable(gl.GL_NORMALIZE)
        gl.glColor4f(color[0], color[1], color[2], color[3])
        glu.gluSphere(self._quadric, 1.0, 32, 16)
        gl.glPopAttrib()
        gl.glPopMatrix()
