"""
main window of the tet viewer, subclasses QMainWindow

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
import pathlib
import csv
import numpy as np

import PyQt5.QtWidgets as qw
import PyQt5.QtCore as qc

import tetmeshtools.meshtools.tetgenread as tr
import tetmeshtools.meshtools.tetgenstructs as ts
import tetmeshtools.meshtools.ffeavolfilereader as fr
import tetmeshtools.tetprops as tp
import tetmeshtools.meshtools.tetmodel as tm
import tetmeshtools.meshtools.trisurface as tris
import tetmeshtools.meshtools.tetmesh as tmes

from tetmeshtools.app_tgv.gui.Ui_tetgenviewermain import Ui_TetgenViewerMain

class TetgenViewerMain(qw.QMainWindow, Ui_TetgenViewerMain):
    """main window of the tet viewer, subclasses QMainWindow"""

    ## the title for the window
    window_title = "TetgenViewer"

    def __init__(self,
                 parent=None,
                 config_args=None):
        """
        the object initalization function
            Args:
                parent (QObject): the parent QObject for this window
                config_args (argparse.Namespace): command line initial setup arguments
        """
        super().__init__(parent)
        self.setupUi(self)

        ## model
        self._model = None

        ## connect up signals
        self._tetViewer.reset_input.connect(self.reset_view)

        ## current source directory
        self._current_source = pathlib.Path.home()

        if config_args is None or config_args.input is None:
            return

        if config_args.input.suffix == ".vol":
            self.load_ffea_file(config_args.input)
        else:
            self.load_tetgen_files(config_args.input)

    def load_ffea_file(self, file_path):
        """
        load a ffea .vol file
        Args:
            file_path (pathlib.Path): file path
        """
        if not file_path.exists():
            qw.QMessageBox.warning(
                self,
                self.tr('TetgenViewer'),
                self.tr(f"File {file_path} doesn't exist!"))
            return

        try:
            points, surface, volume = fr.read_file(file_path)

            nodes = {}
            for index, point in enumerate(points):
                index += 1
                nodes[index] = ts.NodePoint(index, point[0], point[1], point[2])

            faces = {}
            for index, face in enumerate(surface):
                index += 1
                faces[index] = ts.Face(index, face[0], face[1], face[2], -1)

            tets = {}
            for index, tet in enumerate(volume):
                index += 1
                tets[index] = ts.Tetrahedron4(index, tet[0], tet[1], tet[2], tet[3], None)

            tet_props = tp.get_tet_props(nodes, tets)

            self._model = tm.TetModel(tris.TriSurface(nodes, faces),
                                      tmes.TetMesh(nodes, tets))

            self.list_tets(tet_props)
            self.display_total_volume(tet_props)
            self.display_total_surface_area()
            self.reset_view()
            self._current_source = file_path.parent
            self._tetViewer.reset_view()
            self._tetViewer.set_model(self._model)
            self.set_title(file_path)

        except ValueError as error:
            qw.QMessageBox.warning(self, "Tetgen viewer", error)
            return

    def load_tetgen_files(self, root_name):
        """
        load tetgen files
        Args:
            root_name (pathlib.Path): the path and root name of files
        """
        tet_files = []
        tet_files.append(root_name.with_suffix(".1.node"))
        tet_files.append(root_name.with_suffix(".1.face"))
        tet_files.append(root_name.with_suffix(".1.ele"))

        for file in tet_files:
            if not file.exists():
                qw.QMessageBox.warning(
                    self,
                    self.tr('TetgenViewer'),
                    self.tr(f"File {file} doesn't exist!"))
                return

        try:
            _, nodes = tr.read_node_file(tet_files[0])
            _, faces = tr.read_face_file(tet_files[1])
            _, tets  = tr.read_tet_file(tet_files[2])

            tet_props = tp.get_tet_props(nodes, tets)

            self._model = tm.TetModel(tris.TriSurface(nodes, faces),
                                      tmes.TetMesh(nodes, tets))

            self.list_tets(tet_props)
            self.display_total_volume(tet_props)
            self.display_total_surface_area()
            self.reset_view()
            self._current_source = root_name.parent
            self._tetViewer.set_model(self._model)
            self.set_title(root_name)

        except ValueError as error:
            qw.QMessageBox.warning(self, "Tetgen viewer", error)
            return

    def set_title(self, file_path=None):
        """
        add the file name to the window title
        Args:
            file_path (pathlib.Path): the file
        """
        if file_path is None:
            self.setWindowTitle(TetgenViewerMain.window_title)

        self.setWindowTitle(TetgenViewerMain.window_title + ' ' + file_path.name)

    def display_total_volume(self, tet_props):
        """
        display the total volume
        Args:
            tet_props
        """
        total = sum(props[2] for props in tet_props.values())
        self._totalVolLabel.setText(str(round(total, 2)))

    def display_total_surface_area(self):
        """
        find and display the total surface area of the model
        """
        area = self._model.get_surface().surface_area()
        self._totalAreaLabel.setText(str(round(area, 2)))

    def reset_view(self):
        """
        on load reset zero rotations, shift and prespective
        """
        self.reset_sliders()

        old_state = self._surfaceLatticeButton.blockSignals(True)
        self._surfaceLatticeButton.setChecked(False)
        self._surfaceLatticeButton.blockSignals(old_state)

        old_state = self._surfaceButton.blockSignals(True)
        self._surfaceButton.setChecked(False)
        self._surfaceButton.blockSignals(old_state)

        old_state = self._rotXSlider.blockSignals(True)
        self._rotXSlider.setSliderPosition(0)
        self._rotXSlider.blockSignals(old_state)

        old_state = self._rotYSlider.blockSignals(True)
        self._rotYSlider.setSliderPosition(0)
        self._rotYSlider.blockSignals(old_state)

        old_state = self._thicknessBox.blockSignals(True)
        self._thicknessBox.setValue(1)
        self._thicknessBox.blockSignals(old_state)

        self._tetViewer.reset_view()

    def list_tets(self, tet_props):
        """
        list the tet properties in the table
        Args:
            tet_props (dict(int, [float])): tetgen index => properties of the tet
        """
        headers = ["TetGen Index", "Shortest Side", "Volume", "Surface Area", "Shape Factor"]

        old_state = self._tetsTableWidget.blockSignals(True)

        self._tetsTableWidget.clear()
        self._tetsTableWidget.setHorizontalHeaderLabels(headers)
        self._tetsTableWidget.setRowCount(len(tet_props))
        self._tetsTableWidget.setSortingEnabled(True)
        self._tetsTableWidget.verticalHeader().hide()

        for row, key_value in enumerate(tet_props.items()):
            item = qw.QTableWidgetItem()
            item.setData(qc.Qt.DisplayRole, key_value[0])
            self._tetsTableWidget.setItem(row, 0, item)
            item = qw.QTableWidgetItem()
            item.setData(qc.Qt.DisplayRole, float(round(key_value[1][0], 2)))
            self._tetsTableWidget.setItem(row, 1, item)
            item = qw.QTableWidgetItem()
            item.setData(qc.Qt.DisplayRole, float(round(key_value[1][1], 2)))
            self._tetsTableWidget.setItem(row, 2, item)
            item = qw.QTableWidgetItem()
            item.setData(qc.Qt.DisplayRole, float(round(key_value[1][2], 2)))
            self._tetsTableWidget.setItem(row, 3, item)
            item = qw.QTableWidgetItem()
            tmp = key_value[1][2]/np.power(key_value[1][1], (2.0/3.0))
            item.setData(qc.Qt.DisplayRole, float(round(tmp, 2)))
            self._tetsTableWidget.setItem(row, 4, item)

        self._tetsTableWidget.blockSignals(old_state)

    @qc.pyqtSlot()
    def get_and_load_files(self):
        """
        callback for loading a file
        """
        path = str(pathlib.Path.home())
        if self._current_source is not None:
            path = str(self._current_source)

        file_types = "node (*.1.node);; face (*.1.face);; element (*.1.ele);; ffea vol (*.vol)"
        fname, ftype = qw.QFileDialog.getOpenFileName(self,
                                                  "Enter one tetgen file",
                                                  path,
                                                  file_types)
        if fname is None or fname == '':
            return

        if ftype == "ffea vol (*.vol)":
            self.load_ffea_file(pathlib.Path(fname))
            return

        for suffix in [".1.node", ".1.ele", ".1.face"]:
            if fname.endswith(suffix):
                fname = fname.replace(suffix, "", 1)

        self.load_tetgen_files(pathlib.Path(fname))

    @qc.pyqtSlot()
    def selection_changed(self):
        """
        callback for change of selection in table
        """
        indices = self._tetsTableWidget.selectedIndexes()
        row = indices[0].row()
        index = self._tetsTableWidget.item(row, 0)

        self._tetViewer.display_tet(int(index.text()))

    @qc.pyqtSlot()
    def save_tet_data(self):
        """
        save table of tet data
        """
        path = str(pathlib.Path.home())
        if self._current_source is not None:
            path = str(self._current_source)

        name, _ = qw.QFileDialog.getSaveFileName(self,
                                                 "Enter one tetgen file",
                                                 path,
                                                 "Commer sep (*.csv)")

        if name is None or name == '':
            return

        file_path = pathlib.Path(name)

        columns = self._tetsTableWidget.columnCount()
        headers = []
        for column in range(columns):
            headers.append(self._tetsTableWidget.horizontalHeaderItem(column).text())

        with file_path.open('w', newline='', encoding='utf-8') as file:
            writer = csv.writer(file, delimiter=',')
            writer.writerow(headers)
            for row in range(self._tetsTableWidget.rowCount()):
                items = []
                for column in range(columns):
                    items.append(self._tetsTableWidget.item(row, column).text())
                writer.writerow(items)

        qw.QMessageBox.information(self,
                                   "TeggenView Save",
                                   f"Data written to {name}")

    @qc.pyqtSlot()
    def save_image(self):
        """
        save image data
        """
        path = str(pathlib.Path.home())
        if self._current_source is not None:
            path = str(self._current_source)

        file_types = "Portable Network Graphics (*.png);;Joint Photographic Experts Group (*.jpg)"
        name, _ = qw.QFileDialog.getSaveFileName(self,
                                                 "Select or enter file for image",
                                                 path,
                                                 file_types)

        if name is None or name == '':
            return

        self._tetViewer.grab().save(name)

        qw.QMessageBox.information(self,
                                   "TetgenView Save",
                                   f"Image written to {name}")

    @qc.pyqtSlot()
    def reset_sliders(self):
        """
        reset the slider to zero
        """
        old_state = self._rotXSlider.blockSignals(True)
        self._rotXSlider.setSliderPosition(0)
        self._rotXSlider.blockSignals(old_state)

        old_state = self._rotYSlider.blockSignals(True)
        self._rotYSlider.setSliderPosition(0)
        self._rotYSlider.blockSignals(old_state)

    @qc.pyqtSlot()
    def save_setup(self):
        """
        save the current setup
        """
        if self._model is None:
            qw.QMessageBox.information(self,
                                       "Error",
                                       "You must load a mesh to have a setup")
            return

        path = str(pathlib.Path.home())
        name, _ = qw.QFileDialog.getSaveFileName(self,
                                                 "Enter file",
                                                 path,
                                                 "json (*.json)")

        if name is None or name == '':
            return

        file_path = pathlib.Path(name)
        self._tetViewer.save_setup(file_path)

    @qc.pyqtSlot()
    def load_setup(self):
        """
        load a setup file
        """
        if self._model is None:
            qw.QMessageBox.information(self,
                                       "Error",
                                       "You must load a mesh to have a setup")
            return

        path = str(pathlib.Path.home())
        name, _ = qw.QFileDialog.getOpenFileName(self,
                                                 "Enter file",
                                                 path,
                                                 "json (*.json)")

        if name is None or name == '':
            return

        file_path = pathlib.Path(name)
        self._tetViewer.load_setup(file_path)

    @qc.pyqtSlot(qw.QAbstractButton, bool)
    def background_change(self, button, flag):
        """
        callback for selection of a new background
        """
        if not flag:
            return

        self._tetViewer.change_background(button.text())

    @qc.pyqtSlot()
    def view_perspective(self):
        """
        responde to request for prespective view
        """
        self._tetViewer.set_view("Perspective")

    @qc.pyqtSlot()
    def view_orthogonal(self):
        """
        responde to request for orthogonal view
        """
        self._tetViewer.set_view("Orthogonal")

    @qc.pyqtSlot()
    def background_black(self):
        """
        responde to request for black background
        """
        self._tetViewer.change_background("Black")

    @qc.pyqtSlot()
    def background_white(self):
        """
        responde to request for white background
        """
        self._tetViewer.change_background("White")

    @qc.pyqtSlot()
    def background_gray(self):
        """
        responde to request for gray background
        """
        self._tetViewer.change_background("Gray")
