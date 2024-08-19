import unittest
import pathlib
import zipfile
from io import TextIOWrapper

import tetmeshtools.meshtools.tetgenread as tr
import tetmeshtools.meshtools.tetgenwrite as tw

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

class TestWriteTetgen(unittest.TestCase):
    """
    tests of Tetgen file write
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
        self.cleanup = list()
    
    def tearDown(self):
        """
        clean up
        """
        for path in self.cleanup:
            path.unlink()

    def test_write_nodes(self):
        """
        test writing the nodes files
        """
        # Load initial values
        with self._archive.open('sphere10.1.node', 'r') as file:
            _, nodes = tr.read_node_text(TextIOWrapper(file))
        
        # Save copy of values
        tw.write_tetgen_nodes(pathlib.Path("tetmeshtools\\tests\\data\\__nodes_temp"), nodes)

        copy_path = pathlib.Path("tetmeshtools\\tests\\data\\__nodes_temp.1.node")

        # Read copy of values
        _, nodes_copy = tr.read_node_file(copy_path)

        # Add file to cleanup.
        self.cleanup.append(copy_path)

        # Verify that they are identical.
        self.assertEqual(nodes, nodes_copy)

    def test_write_faces(self):
        """
        test writing the faces files
        """
        # Load initial values
        with self._archive.open('sphere10.1.face', 'r') as file:
            _, faces = tr.read_face_text(TextIOWrapper(file))
        
        # Save copy of values
        tw.write_tetgen_faces(pathlib.Path("tetmeshtools\\tests\\data\\__faces_temp"), faces)

        copy_path = pathlib.Path("tetmeshtools\\tests\\data\\__faces_temp.1.face")

        # Read copy of values
        _, faces_copy = tr.read_face_file(copy_path)

        # Add file to cleanup.
        self.cleanup.append(copy_path)

        # Verify that they are identical.
        self.assertEqual(faces, faces_copy)

    def test_write_tets(self):
        """
        test writing the element files
        """
        # Load initial values
        with self._archive.open('sphere10.1.ele', 'r') as file:
            _, tets = tr.read_tet_text(TextIOWrapper(file))
        
        # Save copy of values
        tw.write_tetgen_elements(pathlib.Path("tetmeshtools\\tests\\data\\__tets_temp"), tets)

        copy_path = pathlib.Path("tetmeshtools\\tests\\data\\__tets_temp.1.ele")

        # Read copy of values
        _, tets_copy = tr.read_tet_file(copy_path)

        # Add file to cleanup.
        self.cleanup.append(copy_path)

        # Verify that they are identical.
        self.assertEqual(tets, tets_copy)
