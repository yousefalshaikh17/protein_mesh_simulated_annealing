"""
unit tests for the tetmeshviewr

-----------------------------------------------

Licensed under the GNU General Public License, Version 3.0 (the "License"); you
may not use this file except in compliance with the License. You may obtain a
copy of the License at <https://www.gnu.org/licenses/gpl-3.0.html>.

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

This work was funded by Joanna Leng's EPSRC funded RSE Fellowship (EP/R025819/1)

@copyright 2024
@author: j.h.pickering@leeds.ac.uk and j.leng@leeds.ac.uk
"""
import unittest
import sys

import PyQt6.QtWidgets as qw
from PyQt6.QtCore import Qt
from PyQt6 import QtTest

# set up linting
# pylint: disable = import-error
# pylint: disable = protected-access
# pylint: disable = c-extension-no-member

from tetmeshtools.app_tgv.gui.tetgenviewermain import TetgenViewerMain

app = qw.QApplication(sys.argv)

class TestGuiControls(unittest.TestCase):
    """
    test the video control widget
    """

    def setUp(self) -> None:
        self.app = app.instance()

    def tearDown(self) -> None:
        self.app.closeAllWindows()

    def test_initial_state(self):
        """
        test initialized ok
        """
        main_win = TetgenViewerMain()
        self.assertEqual(main_win._rotXSlider.value(), 0,
                         "X-rotation Slider not initialized correctly")
        self.assertEqual(main_win._rotYSlider.value(), 0,
                         "Y-rotation Slider not initialized correctly")

    def test_sliders_limits(self):
        """
        test motion of sliders.
        """
        main_win = TetgenViewerMain()
        self.assertEqual(main_win._rotXSlider.maximum(), 360,
                         "X-rotation Slider maximum not 360")
        self.assertEqual(main_win._rotXSlider.minimum(), 0,
                         "X-rotation Slider minimum not 0")
        self.assertEqual(main_win._rotYSlider.maximum(), 360,
                         "Y-rotation Slider maximum not 360")
        self.assertEqual(main_win._rotYSlider.minimum(), 0,
                         "Y-rotation Slider minimum not 0")

    def test_check_boxes_initial_states(self):
        """
        test the checkboxes.
        """
        main_win = TetgenViewerMain()
        self.assertTrue(main_win._showTetBox.isChecked(),
                        "Show tet box bad initial state")
        self.assertFalse(main_win._surfaceLatticeButton.isChecked(),
                         'Show lattice is initially checked')

    def test_showTetBox(self):
        """
        test the checkboxes.
        """
        main_win = TetgenViewerMain()

        spy = QtTest.QSignalSpy(main_win._showTetBox.stateChanged)
        QtTest.QTest.mouseClick(main_win._showTetBox, Qt.MouseButton.LeftButton)
        self.assertEqual(len(spy), 1, "_showTetBox: the wrong number of signals emitted")
        self.assertEqual(spy[0][0], Qt.CheckState.Unchecked.value, "_showTetBox wrong check state")
        self.assertFalse(main_win._tetViewer._state._display_current_tet, "viewer not changing _dispaly_current_tet")

    def test_surfaceButton(self):
        """
        test the checkboxes.
        """
        main_win = TetgenViewerMain()

        spy = QtTest.QSignalSpy(main_win._surfaceButton.stateChanged)
        QtTest.QTest.mouseClick(main_win._surfaceButton, Qt.MouseButton.LeftButton)
        self.assertEqual(len(spy), 1, "_surfaceButton: the wrong number of signals emitted")
        self.assertEqual(spy[0][0], Qt.CheckState.Checked.value, "_surfaceButton wrong check state")
        self.assertTrue(main_win._tetViewer._state._show_faces, "viewer not changing _show_faces")

    def test_surfaceLatticeButton(self):
        """
        test the checkboxes.
        """
        main_win = TetgenViewerMain()

        spy = QtTest.QSignalSpy(main_win._surfaceLatticeButton.stateChanged)
        QtTest.QTest.mouseClick(main_win._surfaceLatticeButton, Qt.MouseButton.LeftButton)
        self.assertEqual(len(spy), 1, "_surfaceLatticeButton: the wrong number of signals emitted")
        self.assertEqual(spy[0][0], Qt.CheckState.Checked.value, "_surfaceLatticeButton wrong check state")
        self.assertTrue(main_win._tetViewer._state._show_lattice, "viewer not changing _show_lattice")
