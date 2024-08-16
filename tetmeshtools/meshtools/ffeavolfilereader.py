"""
Reader for a mesh stored in the Fluctuating Finite Element Analysis (FFEA) input
file format (https://ffea.bitbucket.io/).

--------------------

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
# set up linting conditions
# pylint: disable = import-error
import tetmeshtools.meshtools.tetgenstructs as ts
import tetmeshtools.meshtools.trisurface as trs
import tetmeshtools.meshtools.tetmesh as tm
import tetmeshtools.meshtools.tetmodel as tmod

def make_model_from_ffea(mesh_file):
    """
    read a ffea .vol file and construct a mesh model
    Args:
        mesh_file (pathlib.Path): the source file
    Returns
        (TetModel)
    Throws:
        ValueError if bad file
    """
    points, surface, volume = read_file(mesh_file)

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

    surface = trs.TriSurface(nodes, faces)
    mesh    = tm.TetMesh(nodes, tets)

    return tmod.TetModel(surface, mesh, mesh_file)

def get_data_start_end_indices(lines, keyword):
    """
    get the start and end indices for a block beginning with keyword
    Args:
        lines (str): source text, newlines stripped
        keywork (str): text marker for start of block
    Returns:
        tuple(int, int): the start and end indices of data
    """
    key_line = lines.index(keyword)
    number_entries = int(lines[key_line+1])
    start = key_line+2

    return (start, start+number_entries)

def get_points(lines):
    """
    extract the points from the file contents
    Args:
        lines ([str]): text file lines with endlines removed
    Returns
        [[float, float, float]]: the points
    """
    start, end = get_data_start_end_indices(lines, "points")
    points = []
    for index in range(start, end):
        parts = lines[index].split()
        points.append([float(x) for x in parts])

    return points

def get_surface(lines):
    """
    extract the triangles from the file contents
    Args:
        lines ([str]): text file lines with endlines removed
    Returns
        [[int, int, int]]: the points
    """
    start, end = get_data_start_end_indices(lines, "surfaceelements")
    triangles = []
    for index in range(start, end):
        parts = lines[index].split()
        triangles.append([int(parts[5]), int(parts[6]), int(parts[7])])

    return triangles

def get_volume(lines):
    """
    extract the tets from the file contents
    Args:
        lines ([str]): text file lines with endlines removed
    Returns
        [[int, int, int, int]]: the points
    """
    start, end = get_data_start_end_indices(lines, "volumeelements")
    tets = []
    for index in range(start, end):
        parts = lines[index].split()
        tets.append([int(parts[2]), int(parts[3]), int(parts[4]), int(parts[5])])

    return tets

def read_file(source_file):
    """
    read ffea source file
    Args:
        source_file (pathlib.Path): the file
    Returns
        [[float, float, float]]: vertices
        [[int, int, int]]: surface triangles
        [[int, int, int, int]]: volume tets
    """
    if not source_file.exists():
        raise ValueError(f"File {source_file} does not exist.")

    with source_file.open(mode='r') as in_file:
        lines = in_file.readlines()

    lines = [x.strip() for x in lines]

    points = get_points(lines)
    surface = get_surface(lines)
    volume = get_volume(lines)

    return points, surface, volume
