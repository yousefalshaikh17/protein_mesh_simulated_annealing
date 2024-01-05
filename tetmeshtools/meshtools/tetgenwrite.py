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

def write_tetgen_elements(output_file, tets, comment=""):
    """
    write tetgen elements file .ele
    First line: <# of tetrahedra> <nodes per tetrahedron> <# of attributes>
    Remaining lines list of # of tetrahedra:
    <tetrahedron #> <node> <node> <node> <node> ... [attributes]
    Args:
        output_file (pathlib.Path): root name of file, will have .node added
        tets ({index => Tetrahedron4}):
        comment (str): any user comment
    """
    with open(output_file.with_suffix(".1.ele"), "w", encoding="utf-8") as ele:
        ele.write(f"{len(tets)} 4 0\n")
        for index in tets:
            tet = tets[index]
            ele.write(f"{index} {tet.vert0} {tet.vert1} {tet.vert2} {tet.vert3}\n")
        ele.write(f"# {comment}")

def write_tetgen_nodes(output_file, points, comment=""):
    """
    write tetgen .node file
    First line: <num points> <dimension (3)> <num attributes> <num boundary markers (0 or 1)>
    Remaining lines list # of points:
    <point #> <x> <y> <z>
    Args:
        output_file (pathlib.Path): root name of file, will have .node added
        points ({index => NodePoint}):
        comment (str): any user comment
    """
    with open(output_file.with_suffix(".1.node"), "w", encoding="utf-8") as node:
        node.write(f"{len(points)} 3 0 0\n")
        for index in points:
            point = points[index]
            node.write(f"{index} {point.x} {point.y} {point.z}\n")
        node.write(f"# {comment}")

def write_tetgen_faces(output_file, faces, comment=""):
    """
    write .face file
    First line: <# of faces> <boundary marker (0 or 1)>
    Remaining lines list of # of faces:
    <face #> <node> <node> <node> [boundary marker]
    Args:
        output_file (pathlib.Path): root name of file, will have .face added
        faces ({index => Face}): connectivity of face polygons
        comment (str): any user comment
    """
    with open(output_file.with_suffix(".1.face"), "w", encoding="utf-8") as face_file:
        face_file.write(f"{len(faces)} 1\n")
        for index in faces:
            face = faces[index]
            face_file.write(f"{index} {face.vert0} {face.vert1} {face.vert2} -1\n")
        face_file.write(f"# {comment}")
