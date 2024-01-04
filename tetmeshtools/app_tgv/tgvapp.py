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

import PyQt5.QtWidgets as qw

from tetmeshtools.app_tgv.gui.tetgenviewermain import TetgenViewerMain

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
