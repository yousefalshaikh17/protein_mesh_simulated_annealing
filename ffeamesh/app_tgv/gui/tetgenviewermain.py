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
import pathlib
import csv
import numpy as np

import PyQt5.QtWidgets as qw
import PyQt5.QtCore as qc

import ffeamesh.tetmeshtools.tetgenread as tr
import ffeamesh.tetmeshtools.tetgenstructs as ts
import ffeamesh.tetmeshtools.ffeavolfilereader as fr
import ffeamesh.tetprops as tp
import ffeamesh.tetmeshtools.tetmodel as tm
import ffeamesh.tetmeshtools.trisurface as tris
import ffeamesh.tetmeshtools.tetmesh as tmes

from ffeamesh.app_tgv.gui.Ui_tetgenviewermain import Ui_TetgenViewerMain

class TetgenViewerMain(qw.QMainWindow, Ui_TetgenViewerMain):
    """the viewers main window"""

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

        if config_args.input is None:
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
            self.reset_view()
            self._current_source = file_path.parent
            self._tetViewer.reset_view()
            self._tetViewer.set_model(self._model)

        except ValueError as error:
            qw.QMessageBox.warning(self, "Tetgen viewer", error)
            return

    def load_tetgen_files(self, root_name):
        """
        load tetgen files
        Args:
            root_name (pathlib.Path): the path and root name of files
        """
        node_file = root_name.with_suffix(".1.node")
        face_file = root_name.with_suffix(".1.face")
        tets_file = root_name.with_suffix(".1.ele")

        try:
            _, nodes = tr.read_node_file(node_file)
            _, faces = tr.read_face_file(face_file)
            _, tets  = tr.read_tet_file(tets_file)

            tet_props = tp.get_tet_props(nodes, tets)

            self._model = tm.TetModel(tris.TriSurface(nodes, faces),
                                      tmes.TetMesh(nodes, tets))

            self.list_tets(tet_props)
            self.reset_view()
            self._current_source = root_name.parent
            self._tetViewer.set_model(self._model)

        except ValueError as error:
            qw.QMessageBox.warning(self, "Tetgen viewer", error)
            return

    def reset_view(self):
        """
        on load reset zero rotations, shift and prespective
        """
        self.reset_sliders()

        old_state = self._viewGroup.blockSignals(True)
        self._perspectiveButton.setChecked(True)
        self._viewGroup.blockSignals(old_state)

        old_state = self._backgroundGroup.blockSignals(True)
        self._whiteButton.setChecked(True)
        self._backgroundGroup.blockSignals(old_state)

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
                fname = fname.removesuffix(suffix)

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

    @qc.pyqtSlot(bool)
    def show_surface(self, flag):
        """
        callback for the show model surface widget
        Args:
            flag (bool): the new state
        """
        if self._model is None:
            return

        self._tetViewer.show_faces(flag)

    @qc.pyqtSlot()
    def save_setup(self):
        """
        save the current setup
        """
        # test if there is a set up
        # request outfile
        print("Save")

    @qc.pyqtSlot()
    def load_setup(self):
        """
        load a setup file
        """
        # request infile
        print("load")

    @qc.pyqtSlot(bool)
    def show_surface_lattice(self, flag):
        """
        callback for the show model surface edges widget
        Args:
            flag (bool): the new state
        """
        if self._model is None:
            return

        self._tetViewer.show_surface_lattice(flag)

    @qc.pyqtSlot(qw.QAbstractButton, bool)
    def background_change(self, button, flag):
        """
        callback for the motion of a mouse with a button held
        """
        if not flag:
            return

        self._tetViewer.change_background(button.text())

    @qc.pyqtSlot(qw.QAbstractButton, bool)
    def view_change(self, button, flag):
        """
        callback for the motion of a mouse with a button held
        """
        if not flag:
            return

        self._tetViewer.set_view(button.text())
