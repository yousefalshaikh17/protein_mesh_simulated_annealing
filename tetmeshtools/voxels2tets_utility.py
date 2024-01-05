"""
You should have received a copy of the GNU General Public License.
If not, see <http://www.gnu.org/licenses/>.

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

This work was funded by Joanna Leng's EPSRC funded RSE Fellowship (EP/R025819/1)

@copyright 2020
@author: j.h.pickering@leeds.ac.uk and j.leng@leeds.ac.uk

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
"""
# set up linting
# pylint: disable = import-error
# the following avoids problem with vtk's c++
# pylint: disable = no-member
import datetime
import getpass
import enum
import numpy as np

import vtk
from vtkmodules.util import numpy_support
from tetmeshtools import vtk_write
import tetmeshtools.coord_utility as cu
from tetmeshtools.tetgen_write import write_tetgen_output
import tetmeshtools.vtk_utility as vtk_u

class PruneLevel(enum.IntEnum):
    """
    number of verticies allowed below isovlue
    """
    FOUR  = 4
    THREE = 3
    TWO   = 2
    ONE   = 1

def make_fractional_to_cartesian_conversion_function(mrc):
    """
    make_fractional_to_cartesian_conversion_function(mrc)
    make a functor that will convert a point defined by fractional
    coordinates in the unit cell to cartesian coordinates
    assumes angles are all 90 degrees
    Args:
        mrc (mrcfile): the mrc file
    Returns
        function: fractional (i, j, k) (float) => cartesian (x, y, z) (float)
    """
    # get cell sizes
    x_res = mrc.header.cella.x/np.float32(mrc.header.mx)
    y_res = mrc.header.cella.y/np.float32(mrc.header.my)
    z_res = mrc.header.cella.z/np.float32(mrc.header.mz)

    # get coordinate origine
    x_trans = np.float32(mrc.header.origin['x'])
    y_trans = np.float32(mrc.header.origin['y'])
    z_trans = np.float32(mrc.header.origin['z'])

    def fractional_to_cartesian_coordinates(x_frac, y_frac, z_frac):
        """
        convert integer indices to 3D coordinates.
        Args:
            x_frac (float): the fractional x coordinate
            y_frac (float): the fractional y coordinate
            z_frac (float): the fractional z coordinate
        Returns:
            (CoordTransform): container for fractional and cartesian coordinates
        """
        x_cart = (x_frac * x_res) + x_trans
        y_cart = (y_frac * y_res) + y_trans
        z_cart = (z_frac * z_res) + z_trans

        return cu.CoordTransform(cu.Coordinate(x_frac, y_frac, z_frac),
                                 cu.Coordinate(x_cart, y_cart, z_cart))

    return fractional_to_cartesian_coordinates

def cube_6_tet_indices():
    """
    return a list of lists for the constuction of 6 tets from the 8 vertices a cube
    Args:
        None
    Returns
        lits(list) : five lists of four vertex indices representing the tets
    """
    return [[0, 6, 5, 1],
            [0, 6, 1, 2],
            [0, 6, 2, 3],
            [0, 6, 3, 7],
            [0, 6, 7, 4],
            [0, 6, 4, 5]]

def even_cube_tet_indices():
    """
    return a list of lists for the constuction of 5 tets from the 8 vertices of an even cube
    Args:
        None
    Returns
        lits(list) : five lists of four vertex indices representing the tets
    """
    return [[0, 4, 5, 7],
            [0, 1, 2, 5],
            [2, 5, 6, 7],
            [0, 2, 3, 7],
            [0, 2, 5, 7]]

def even_cube_tets(cube):
    """
    convert a list of the eight vertices of a an even cube
    into a list of lists representing tets
    Args:
        cube (list): the eight vertices of a cube
    Returns
        lits(list) : five lists of four vertices representing the tets
    """
    indices_array = even_cube_tet_indices()
    tet_list = []
    for indices in indices_array:
        tet_list.append(np.array([cube[indices[0]],
                                  cube[indices[1]],
                                  cube[indices[2]],
                                  cube[indices[3]]]))

    return tet_list

def odd_cube_tet_indices():
    """
    return a list of lists for the constuction of 5 tets from the 8 vertices of an odd cube
    Args:
        None
    Returns
        lits(list) : five lists of four vertex indices representing the tets
    """
    return [[0, 1, 3, 4],
            [1, 4, 5, 6],
            [1, 2, 3, 6],
            [3, 4, 6, 7],
            [1, 3, 4, 6]]

def odd_cube_tets(cube):
    """
    convert a list of the eight vertices of a an odd cube
    into a list of lists representing tets
    Args:
        cube (list): the eight vertices of a cube
    Returns
        lits(list) : five lists of four vertices representing the tets
    """
    indices_array = odd_cube_tet_indices()
    tet_list = []
    for indices in indices_array:
        tet_list.append(np.array([cube[indices[0]],
                                  cube[indices[1]],
                                  cube[indices[2]],
                                  cube[indices[3]]]))

    return tet_list

def is_odd(x_index, y_index, z_index):
    '''
    Logic to decide if a voxel is an odd or and even.
    Args:
        x_index (int) Index for x value in values array.
        y_index (int) Index for y value in values array.
        z_index (int) Index for z value in values array.
    Returns:
        (bool)  True if voxel in an odd position, else odd
    '''
    flag = None
    y_parity =  y_index % 2 == 0
    x_parity =  x_index % 2 == 0

    # Logic for alternating tet division 0 (even) or 1 (odd) - to identify the odd and even voxels
    if z_index % 2 == 0:
        if y_parity:
            flag = not x_parity
        else:
            flag = x_parity
    else:
        if y_parity:
            flag = x_parity
        else:
            flag = not x_parity

    return flag

def create_cube_coords(voxel, frac_to_cart):
    '''
    Caluculates the next 8 coords for the next voxel that has been
    previously thresholded (logic in loop that calls this one).
    Args:
        voxel (Coordinate): array indices of voxel
        frac_to_cart (float*3=>CoordTransform): fractional to cartesian conversin function
    Returns:
        None
    '''
    coords = []
    # Calculate the cordinates of each vertex of the voxel
    coords.append(frac_to_cart((voxel.x-0.5), (voxel.y-0.5), (voxel.z-0.5)))
    coords.append(frac_to_cart((voxel.x+0.5), (voxel.y-0.5), (voxel.z-0.5)))
    coords.append(frac_to_cart((voxel.x+0.5), (voxel.y+0.5), (voxel.z-0.5)))
    coords.append(frac_to_cart((voxel.x-0.5), (voxel.y+0.5), (voxel.z-0.5)))
    coords.append(frac_to_cart((voxel.x-0.5), (voxel.y-0.5), (voxel.z+0.5)))
    coords.append(frac_to_cart((voxel.x+0.5), (voxel.y-0.5), (voxel.z+0.5)))
    coords.append(frac_to_cart((voxel.x+0.5), (voxel.y+0.5), (voxel.z+0.5)))
    coords.append(frac_to_cart((voxel.x-0.5), (voxel.y+0.5), (voxel.z+0.5)))

    return coords

def make_vertex_values(x_index, y_index, z_index, mrc):
    """
    make an arry of eight values by averaging at the vertices
    Args:
        x_index (int)       Indecies to the x coord in the coords array.
        y_index (int)       Indecies to the y coord in the coords array.
        z_index (int)       Indecies to the z coord in the coords array.
        mrc (mrcfile.mmap)  image source file
    Returns:
        (float list): eight values
    """
    ctr_value = mrc.data[z_index, y_index, x_index]
    values = []

    # compute the eight averages
    values.append((ctr_value + mrc.data[z_index-1, y_index-1, x_index-1])/2.0)
    values.append((ctr_value + mrc.data[z_index-1, y_index-1, x_index+1])/2.0)
    values.append((ctr_value + mrc.data[z_index-1, y_index+1, x_index+1])/2.0)
    values.append((ctr_value + mrc.data[z_index-1, y_index+1, x_index-1])/2.0)
    values.append((ctr_value + mrc.data[z_index+1, y_index-1, x_index-1])/2.0)
    values.append((ctr_value + mrc.data[z_index+1, y_index-1, x_index+1])/2.0)
    values.append((ctr_value + mrc.data[z_index+1, y_index+1, x_index+1])/2.0)
    values.append((ctr_value + mrc.data[z_index+1, y_index+1, x_index+1])/2.0)

    return values

def write_tets_to_files(points_list, tets_connectivity, output_file, tetgen_out, vtk_out):
    """
    Sets it up so data goes to the tetgen format writer and the vtk writer correctly.
    Args:
        points_list (float*3 list): 3d coordinates of the points forming the tets
        tets_connectivity (int*4 list): for each tet the indices of its vertices in the points list
        output_file (pathlib.Path): name stem for ouput files
        tetgen_out (bool): if true write tegen fromat files
        vtk_out (bool): if true write a vtk file
    Returns:
        None
    """
    # convert coords to np.array
    points_np = np.array(points_list)

    # connectivity in vtk
    cells_con = vtk_u.make_vtk_tet_connectivity(tets_connectivity)

    # make the grid (vtk scene)
    vtk_pts = vtk.vtkPoints()

    # work on deep copy
    vtk_pts.SetData(numpy_support.numpy_to_vtk(points_np, deep=True))

    vtk_grid = vtk.vtkUnstructuredGrid() #create unstructured grid
    vtk_grid.SetPoints(vtk_pts) #assign points to grid
    vtk_grid.SetCells(vtk.VTK_TETRA, cells_con) #assign tet cells to grid

    #write vtk file
    if vtk_out:
        vtk_write.vtk_output(vtk_grid, output_file)

    # write tetgen file for tetgen input
    if tetgen_out:
        tetgen_output(vtk_grid, points_np, tets_connectivity, output_file)

def tetgen_output(vtk_grid, points, tets_connectivity, output_file):
    """
    Constructs the faces and calls the tetgen writer to output file.
    Args:
        vtk_grid (vtk.vtkUnstructuredGrid): vtk scene
        points (float np.ndarray): duplicate free list of vertices
        tets_connectivity (int*4 list): for each tet the indices of its vertices in the points list
        output_file (pathlib.Path): name stem of output files
    Returns:
        None
    """
    surf_points, cell_data, cell_count = vtk_u.get_vtk_surface(vtk_grid)

    # make a connectivity into the tets points
    original_ids = np.zeros((len(surf_points),), dtype='int16')
    for pos, point in enumerate(surf_points):
        # index of surface point in the surface faces's points array
        original_ids[pos] = np.where((points==point).all(axis=1))[0]

    # This holds true if all polys are of the same kind, e.g. triangles.
    assert cell_data.GetNumberOfValues()%cell_count==0

    # reshape the cells array to match tet gen output standard
    col_count = cell_data.GetNumberOfValues()//cell_count
    numpy_cells = np.array(cell_data)
    faces = numpy_cells.reshape((-1, col_count))

    #write to tetgen .ele, .node, .face
    date = datetime.datetime.now().strftime("%x")
    write_tetgen_output(output_file,
                      tets_connectivity,
                      points,
                      faces,
                      original_ids,
                      f'# created by {getpass.getuser()} on {date}')
