import unittest

import PyQt5.QtWidgets as qw
import PyQt5.Qt as qt
from PyQt5.QtTest import QTest, QSignalSpy

from tetmeshtools.app_tgv.gui.tetgenviewermain import TetgenViewerMain

class TestGuiControls(unittest.TestCase):
    """
    test the video control widget
    """

    def setUp(self):
        """
        build a full test class
        """
        ## the QApplication
        self.app = qw.QApplication([])

        ## the widget
        self._main = TetgenViewerMain()

    def tearDown(self):
        """
        remove
        """
        del self._main

    def test_initial_state(self):
        """
        test initialized ok
        """
        self.assertEqual(self._main._rotXSlider.value(), 0, "X-rotation Slider not initialized correctly")
        self.assertEqual(self._main._rotYSlider.value(), 0, "Y-rotation Slider not initialized correctly")

    def test_sliders_limits(self):
        """
        test motion of sliders.
        """
        self.assertEqual(self._main._rotXSlider.maximum(), 360, "X-rotation Slider maximum not 360")
        self.assertEqual(self._main._rotXSlider.minimum(), 0, "X-rotation Slider minimum not 0")
        self.assertEqual(self._main._rotYSlider.maximum(), 360, "Y-rotation Slider maximum not 360")
        self.assertEqual(self._main._rotYSlider.minimum(), 0, "Y-rotation Slider minimum not 0")

    # def test_check_boxes_initial_states(self):
    #     """
    #     test the checkboxes.
    #     """
    #     message = "Show tet box bad initial state"
    #     #self.assertTrue(self._main._tetViewer._showTetBox.isChecked(), message)
    #     self.assertFalse(self._main._surfaceLatticeButton.isChecked(), 'Show lattice is initially checked')

    def test_check_boxes(self):
        """
        test the checkboxes.
        """
        old = self._main._tetViewer._show_faces
        self._main._surfaceButton.setChecked(True)
        QTest.mouseClick(self._main._surfaceButton, qt.Qt.LeftButton)
        message = "Show faces button not changing state"
        self.assertEqual(self._main._tetViewer._show_faces, old, message)
