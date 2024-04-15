"""
Licensed under the GNU General Public License, Version 3.0 (the "License"); you
may not use this file except in compliance with the License. You may obtain a
copy of the License at <https://www.gnu.org/licenses/gpl-3.0.html>.

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

This work was funded by Joanna Leng's EPSRC funded RSE Fellowship (EP/R025819/1)

@copyright 2023
@author: j.h.pickering@leeds.ac.uk and j.leng@leeds.ac.uk
"""
import unittest
import pathlib
import zipfile
from io import TextIOWrapper

import tetmeshtools.meshtools.tetgenread as tr

def test_existance(file_path):
    """
    test if file exists and raise ValueError if not
    Args:
        file_path (pathlib.Path)
    Raises:
        ValueError if does not exist
    """
    message = "Test cannot be run as data file {name} doesn't exist!"
    if not file_path.exists():
        raise ValueError(message.format(name=str(file_path)))

class TestReadTetgen(unittest.TestCase):
    """
    tests of Tetgen file reader
    """

    def setUp(self):
        """
        build a full test class
        """
        path = pathlib.Path(pathlib.Path.cwd())
        path = path.joinpath("tetmeshtools\\tests\\data\\sphere10.zip")
        test_existance(path)

        ## path to zip archive of test data
        self._archive = zipfile.ZipFile(path)

    def test_read_nodes(self):
        """
        test reading the nodes files
        """
        with self._archive.open('sphere10.1.node', 'r') as file:
            meta, nodes = tr.read_node_text(TextIOWrapper(file))

            self.assertIsInstance(nodes, dict)
            self.assertEqual(meta.points, len(nodes))
            self.assertEqual(meta.dimension, 3)

    def test_read_tets(self):
        """
        test reading the elements files
        """
        with self._archive.open('sphere10.1.ele', 'r') as file:
            meta, tets = tr.read_tet_text(TextIOWrapper(file))

            self.assertIsInstance(tets, dict)
            self.assertEqual(meta.tets, len(tets))
            self.assertEqual(meta.nodes, 4)

    def test_read_faces(self):
        """
        test reading the elements files
        """
        with self._archive.open('sphere10.1.ele', 'r') as file:
            meta, faces = tr.read_face_text(TextIOWrapper(file))

            self.assertIsInstance(faces, dict)
            self.assertEqual(meta.faces, len(faces))
