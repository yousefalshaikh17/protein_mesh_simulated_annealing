#!/usr/bin/env python
"""
 make examples of the decomposition of cubes into tetrhedron meshs.

      Vertex Indices:
         7+----------+6
         /|         /|
        / |        / |
      4+----------+5 |
       |  |       |  |         Axes:
       | 3+-------|--+2        z  y
       | /        | /          | /
       |/         |/           |/
      0+----------+1           +----x

---------------------------------

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
# set up linting
# pylint: disable = import-error
# pylint: disable = no-member
# pylint: disable = no-name-in-module
import pathlib
import argparse
import numpy as np
import vtk.util.numpy_support

import tetmeshtools.vtk_utility as vtk_u
from tetmeshtools import vtk_write

def unit_cube_coords():
    """
    vertices of unit cube
    Returns:
        ([[float]]): 3x16 the coordinates of the vertices of two cubes
    """
    return [[0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0], # cube 0 end
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 0.0, 1.0],
            [2.0, 0.0, 1.0],
            [2.0, 1.0, 1.0],
            [1.0, 1.0, 1.0]]

def cube_5tet_indices():
    """
    return a list of lists for the constuction of 5 tets from the 8 vertices of an even cube
    Args:
        None
    Returns
        [[int]] : connectivity of five tets in terms of cube vertices
    """
    return [[0, 4, 5, 7],
            [0, 1, 2, 5],
            [2, 5, 6, 7],
            [0, 2, 3, 7],
            [0, 2, 5, 7]]

def cube_6tet_indices():
    """
    return a list of lists for the constuction of 5 tets from the 8 vertices of an even cube
    Args:
        None
    Returns
        [[int]] : connectivity of six tets in terms of vertices of two cubes
    """
    return [[0, 6, 5, 1],
            [0, 6, 1, 2],
            [0, 6, 2, 3],
            [0, 6, 3, 7],
            [0, 6, 7, 4],
            [0, 6, 4, 5],
            [8, 14, 13, 9],
            [8, 14, 9, 10],
            [8, 14, 10, 11],
            [8, 14, 11, 15],
            [8, 14, 15, 12],
            [8, 14, 12, 13]]

def get_args():
    """
    get the command line options
    Returns
        (argparse.namespace)
    """
    description = ("Make a pair of voxels decomposed into tetrahedra "
                   "and output in VTK format. 5 and 6 tetrahedra decompositions "
                   "are available.")
    parser = argparse.ArgumentParser(description=description)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-5",
                       "--five_tets",
                       action='store_true',
                       help="if set decomp into 5 tets")

    group.add_argument("-6",
                       "--six_tets",
                       action='store_true',
                       help="if set decomp into 6 tets")

    return parser.parse_args()

def main():
    """run the script"""
    args = get_args()

    connectivity = None
    file_name = "cube"
    if args.five_tets:
        connectivity = cube_5tet_indices()
        file_name += f"_{5}tets.vkt"
    else:
        connectivity = cube_6tet_indices()
        file_name += f"_{6}tets.vkt"

    points_np = np.array(unit_cube_coords())
    cells_con = vtk_u.make_vtk_tet_connectivity(connectivity)

    # make the grid (vtk scene)
    vtk_pts = vtk.vtkPoints()
    vtk_pts.SetData(vtk.util.numpy_support.numpy_to_vtk(points_np, deep=True))
    grid = vtk.vtkUnstructuredGrid() #create unstructured grid
    grid.SetPoints(vtk_pts) #assign points to grid
    grid.SetCells(vtk.VTK_TETRA, cells_con) #assign tet cells to grid

    # write vtk file
    vtk_write.vtk_output(grid, pathlib.Path(file_name))

if __name__ == '__main__':
    main()
