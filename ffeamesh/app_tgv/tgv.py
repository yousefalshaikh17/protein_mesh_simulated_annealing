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
import argparse
import pathlib
import sys

import PyQt5.QtCore as qc
import PyQt5.QtWidgets as qw

#from ffeamesh.app_tgv.tgvapp import TGVApp
from ffeamesh.app_tgv.gui.tetgenviewermain import TetgenViewerMain

class TGVApp(qw.QApplication):
    """
    the application that runs the main window
    """

    def __init__(self, args, python_args):
        """
        initalize and run main window
            Args:
                args [string] the command line arguments
                python_args (argparse.Namespace) argumens found by argparse
         """
        super().__init__(args)
        self.setApplicationName("TetgenViewerMain")
        self.setApplicationVersion("0.0.0")
        self.setOrganizationName("School of Computer Science, University of Leeds, Leeds, UK")
        self.setOrganizationDomain("leeds.ac.uk")

        window = TetgenViewerMain(config_args=python_args)
        window.resize(500, 300)
        window.show()

        self.exec_()

def get_args():
    """
    get command line args
    """
    parser = argparse.ArgumentParser("""process MRC files and produces a regular
        tetrahedral volumetric mesh for FFEA using the "marching tet" algorithm.
        This is written out in the tetgen .ele, .face, and .node file format for
        later use in FFEA, and .vtk for mesh analysis.

        Joanna Leng (J.Leng@leeds.ac.uk),
        Jonathan Pickering (J.H.Pickering@leeds.ac.uk)""")

    parser.add_argument("-i",
                        "--input",
                        type=pathlib.Path,
                        required=False,
                        help="input file root name")

    return parser.parse_args()

def main():
    """
    run the main window
    """
    # enable for high resolution
    if hasattr(qc.Qt, 'AA_EnableHighDpiScaling'):
        qc.QCoreApplication.setAttribute(qc.Qt.AA_EnableHighDpiScaling, True)
    if hasattr(qc.Qt, 'AA_UseHighDpiPixmaps'):
        qc.QCoreApplication.setAttribute(qc.Qt.AA_UseHighDpiPixmaps, True)

    TGVApp(sys.argv, get_args())

if __name__ == "__main__":
    main()
