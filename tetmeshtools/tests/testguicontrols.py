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
        x_initial = self._main._rotXSlider.value()
        self.assertEqual(x_initial, 0, "Slider not initialized correctly")
