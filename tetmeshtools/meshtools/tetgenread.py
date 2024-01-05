"""
You should have received a copy of the GNU General Public License.
If not, see <http://www.gnu.org/licenses/>.

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

This work was funded by Joanna Leng's EPSRC funded RSE Fellowship (EP/R025819/1)

@copyright 2023
@author: j.h.pickering@leeds.ac.uk and j.leng@leeds.ac.uk
"""
# set up linting conditions
# pylint: disable = import-error
import csv
import tetmeshtools.meshtools.tetgenread as tr
import tetmeshtools.meshtools.tetgenstructs as ts
import tetmeshtools.meshtools.trisurface as trs
import tetmeshtools.meshtools.tetmesh as tm
import tetmeshtools.meshtools.tetmodel as tmod

def make_model_from_tetgen(file_root):
    """
    read tetgen .1.node, .1.face, .1.ele files and construct a mesh model
    Args:
        mesh_file (pathlib.Path): root name of source files
    Returns
        (TetModel)
    Throws:
        ValueError if bad file
    """
    tetgen_files = make_and_test_tetgen_files(file_root)

    _, nodes = tr.read_node_file(tetgen_files[0])
    _, faces = tr.read_face_file(tetgen_files[1])
    _, tets  = tr.read_tet_file(tetgen_files[2])

    surface = trs.TriSurface(nodes, faces)
    mesh    = tm.TetMesh(nodes, tets)

    return tmod.TetModel(surface, mesh, file_root)

def make_and_test_tetgen_files(root):
    """
    make the names of the tetgen files from the name root
    Args:
        root (pathlib.Path): name root
    Returns
        [pathlib.Path]: <name>.1.ele, <name>.1.node, <name>.1.face
    Throws:
        ValueError: if a file does not exist
    """
    files = []
    files.append(root.with_suffix(".1.node"))
    files.append(root.with_suffix(".1.face"))
    files.append(root.with_suffix(".1.ele"))

    for file in files:
        if not file.exists():
            raise ValueError(f"File {file} does not exist!")

    return files

def make_decomment(comment):
    """
    make a decomment function
    Args:
        comment (str): the comment string
    """
    if comment is None or comment.isspace() or comment == '':
        return None

    def decomment(csvfile):
        """
        make iterator to remove lines starting with the comment symbols
        Args:
            csvfile ()
        Returns:
            generator returning non-comment file rows as (str)
        """
        for row in csvfile:
            if len(row)>0 and not row.strip().startswith(comment):
                row = row[:row.find(comment)].strip()
                if len(row)>0:
                    yield row

    return decomment

def read_node_file(input_file):
    """
    read a tetgen nodes file
    Args:
        input_file (pathlib.Path): source
    Return
        (NodeMetaData): the file's metadata
        ({NodePoint}): the points
    Throws:
        ValueError if problem
    """
    with input_file.open('r', newline='') as file:
        return read_node_text(file)

def read_node_text(text_stream):
    """
    read a text stream of tetgen nodes data
    Args:
        text_stream (io.TextIOWrapper)
    """
    decomment = make_decomment("#")

    reader = csv.reader(decomment(text_stream), delimiter=' ')
    row = next(reader)
    row = [x for x in row if x != '']
    meta_data = ts.NodeMetaData(int(row[0]), int(row[1]), int(row[2]), int(row[3]))

    if meta_data.dimension != 3:
        raise ValueError("File input is not 3D.")

    points = {}
    for row in reader:
        row = [x for x in row if x != '']
        index = int(row[0])
        points[index] = ts.NodePoint(index, float(row[1]),  float(row[2]), float(row[3]))

    if len(points) != meta_data.points:
        req = meta_data.points
        act = len(points)
        er_m = f"File input should have {req} points, but {act} were found!"
        raise ValueError(er_m)

    return meta_data, points

def read_face_file(input_file):
    """
    read a tetgen faces file
    Args:
        input_file (pathlib.Path): source
    Return
        (FaceMetaData): the file's metadata
        ([Face]): the faces
    Throws:
        ValueError if problem
    """
    with input_file.open('r', newline='') as file:
        return read_face_text(file)

def read_face_text(text_stream):
    """
    read a text stream of tetgen face data
    Args:
        text_stream (io.TextIOWrapper)
    """
    decomment = make_decomment("#")

    reader = csv.reader(decomment(text_stream), delimiter=' ')
    row = next(reader)
    row = [x for x in row if x != '']
    meta_data = ts.FaceMetaData(int(row[0]), int(row[1]))

    faces = {}
    for row in reader:
        row = [x for x in row if x != '']
        faces[int(row[0])] = ts.Face(int(row[0]),
                                        int(row[1]),
                                        int(row[2]),
                                        int(row[3]),
                                        int(row[4]))

    if len(faces) != meta_data.faces:
        req = meta_data.faces
        act = len(faces)
        er_m = f"File should have {req} faces, but {act} were found!"
        raise ValueError(er_m)

    return meta_data, faces

def read_tet_file(input_file):
    """
    read a tetgen ele file
    Args:
        input_file (pathlib.Path): source
    Return
        (TetMetaData): the file's metadata
        ({index => Tetrahedron4}): the tets
    Throws:
        ValueError if problem
    """
    with input_file.open('r', newline='') as file:
        return read_tet_text(file)

def read_tet_text(text_stream):
    """
    read a text stream of tetgen element data
    Args:
        text_stream (io.TextIOWrapper)
    """
    decomment = make_decomment("#")

    reader = csv.reader(decomment(text_stream), delimiter=' ')
    row = next(reader)
    row = [x for x in row if x != '']
    meta_data = ts.TetMetaData(int(row[0]), int(row[1]), int(row[2]))

    if meta_data.nodes != 4:
        raise ValueError("This reader can only handle four nodes per tetrahedron.")

    tets = {}
    for row in reader:
        row = [x for x in row if x != '']
        index = int(row[0])
        if meta_data.ra == 0:
            tets[index] = ts.Tetrahedron4(index,
                                          int(row[1]),
                                          int(row[2]),
                                          int(row[3]),
                                          int(row[4]),
                                          None)
        else:
            tets[index] = ts.Tetrahedron4(index,
                                          int(row[1]),
                                          int(row[2]),
                                          int(row[3]),
                                          int(row[4]),
                                          int(row[5]))

    if len(tets) != meta_data.tets:
        req = meta_data.tets
        act = len(tets)
        er_m = f"File should have {req} tets, but {act} were found!"
        raise ValueError(er_m)

    return meta_data, tets
