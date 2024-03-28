import unittest
import sys

import PyQt5.QtWidgets as qw
import PyQt5.Qt as qt
from PyQt5.QtTest import QTest, QSignalSpy

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
        self.assertEqual(main_win._rotXSlider.value(), 0, "X-rotation Slider not initialized correctly")
        self.assertEqual(main_win._rotYSlider.value(), 0, "Y-rotation Slider not initialized correctly")

    def test_sliders_limits(self):
        """
        test motion of sliders.
        """
        main_win = TetgenViewerMain()
        self.assertEqual(main_win._rotXSlider.maximum(), 360, "X-rotation Slider maximum not 360")
        self.assertEqual(main_win._rotXSlider.minimum(), 0, "X-rotation Slider minimum not 0")
        self.assertEqual(main_win._rotYSlider.maximum(), 360, "Y-rotation Slider maximum not 360")
        self.assertEqual(main_win._rotYSlider.minimum(), 0, "Y-rotation Slider minimum not 0")

    def test_check_boxes_initial_states(self):
        """
        test the checkboxes.
        """
        main_win = TetgenViewerMain()
        self.assertTrue(main_win._showTetBox.isChecked(), "Show tet box bad initial state")
        self.assertFalse(main_win._surfaceLatticeButton.isChecked(), 'Show lattice is initially checked')

    def test_check_boxes(self):
        """
        test the checkboxes.
        """
        main_win = TetgenViewerMain()

        spy0 = QSignalSpy(main_win._showTetBox.stateChanged)
        spy1 = QSignalSpy(main_win._surfaceButton.stateChanged)
        spy2 = QSignalSpy(main_win._surfaceLatticeButton.stateChanged)

        QTest.mouseClick(main_win._showTetBox, qt.Qt.LeftButton)
        QTest.mouseClick(main_win._surfaceButton, qt.Qt.LeftButton)
        QTest.mouseClick(main_win._surfaceLatticeButton, qt.Qt.LeftButton)

        self.assertEqual(len(spy0), 1, "_showTetBox: the wrong number of signals emitted")
        self.assertEqual(len(spy1), 1, "_surfaceButton: the wrong number of signals emitted")
        self.assertEqual(len(spy2), 1, "_surfaceLatticeButton: the wrong number of signals emitted")

        self.assertEqual(spy0[0][0], qt.Qt.CheckState.Unchecked, "_showTetBox wrong check state")
        self.assertEqual(spy1[0][0], qt.Qt.CheckState.Checked, "_showTetBox wrong check statee")
        self.assertEqual(spy2[0][0], qt.Qt.CheckState.Checked, "_showTetBox wrong check state")
