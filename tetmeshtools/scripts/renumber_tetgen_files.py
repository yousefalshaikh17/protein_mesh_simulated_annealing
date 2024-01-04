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
# from_start_zero
#
# 'fix' tetgen files with nodes starting at zero
import argparse
import pathlib
import itertools

import tetmeshtools.tetmeshtools.tetgenread as tr
import tetmeshtools.tetmeshtools.tetgenwrite as tw
import tetmeshtools.tetmeshtools.tetgenstructs as ts

def  reindex_tets(tets, map):
    """
    reindex the faces to start at one and ascend in steps of one, also use map to update the node indixes
    """
    new_tets = {}
    index = itertools.count(1)

    for tet in tets.values():
        new_tet = ts.Tetrahedron4(next(index),
                                  map[tet.vert0],
                                  map[tet.vert1],
                                  map[tet.vert2],
                                  map[tet.vert3],
                                  tet.ra)
        new_tets[new_tet.index] = new_tet

    return new_tets

def  reindex_faces(faces, map):
    """
    reindex the faces to start at one and ascend in steps of one, also use map to update the node indixes
    """
    new_faces = {}
    index = itertools.count(1)

    for face in faces.values():
        new_face = ts.Face(next(index),
                           map[face.vert0],
                           map[face.vert1],
                           map[face.vert2],
                           face.bm)
        new_faces[new_face.index] = new_face

    return new_faces

def reindex_nodes(old_nodes):
    """
    reindex nodes to start at one and ascend in steps of one
    Args:
        old_nodes dict(index:NodePoint)
    Returns
        dict(old_index => new_index): map of reindexing
        dict(index => node): corrected nodes
    """
    new_nodes = {}
    index_map={}
    index = itertools.count(1)
    for node in old_nodes.values():
        new_node = ts.NodePoint(next(index), node.x, node.y, node.z)
        new_nodes[new_node.index] = new_node
        index_map[node.index] = new_node.index

    return index_map, new_nodes

def existing_tetgen_file(file_path):
    """
    test file_path is a .1.node file make .node, .face & .ele pathlib.Pahts and they all exist
    Args:
        file_path (str)
    Returns
        [pathlib.Path]: .node, .face, .ele
    Raises
        ValueError
    """
    if not str(file_path).endswith(".1.node"):
        raise ValueError(f"file {file_path} is not a .1.node file")

    paths = []
    paths.append(pathlib.Path(file_path))
    paths.append(pathlib.Path(file_path.replace(".1.node", ".1.face")))
    paths.append(pathlib.Path(file_path.replace(".1.node", ".1.ele")))

    for path in paths:
        if not path.exists():
            raise ValueError(f"file {path} doesn't exist")

    paths.append(pathlib.Path(file_path.replace('.1.node', '')))

    return paths

def get_args():
    """
    get command line
    Returns
        namespace
    """
    description = """
renumber tetgen files so that:
    1. the node indices start at one
    2. the node idices increment by one
This makes then compatible with FFEA
"""
    parser = argparse.ArgumentParser(description)

    parser.add_argument('-i',
                        '--input',
                        required=True,
                        type=existing_tetgen_file,
                        help="name of source node file")

    return parser.parse_args()

def main():
    args = get_args()

    _, nodes = tr.read_node_file(args.input[0])
    _, faces = tr.read_face_file(args.input[1])
    _, tets  = tr.read_tet_file(args.input[2])

    map, new_nodes = reindex_nodes(nodes)
    new_faces = reindex_faces(faces, map)
    new_tets = reindex_tets(tets, map)

    nodes_file = args.input[0].rename(pathlib.Path(
                                        args.input[0].parent,
                                        args.input[0].name.replace(".1.node", ".old.1.node")))
    face_file = args.input[1].rename(pathlib.Path(
                                        args.input[1].parent,
                                        args.input[1].name.replace(".1.face", ".old.1.face")))
    tets_file = args.input[2].rename(pathlib.Path(
                                        args.input[2].parent,
                                        args.input[2].name.replace(".1.ele", ".old.1.ele")))

    tw.write_tetgen_nodes(args.input[3], new_nodes, comment="renumbered")
    tw.write_tetgen_faces(args.input[3], new_faces, comment="renumbered")
    tw.write_tetgen_elements(args.input[3], new_tets, comment="renumbered")

if __name__ == '__main__':
    main()
