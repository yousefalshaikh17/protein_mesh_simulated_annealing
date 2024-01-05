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

import sys
from time import process_time
import mrcfile

import tetmeshtools.voxels2tets_utility as v2t
import tetmeshtools.mrclininterpolate as mi
from tetmeshtools.mrc_utility import voxel_size
from tetmeshtools.tetprops import tet_volume
from tetmeshtools.meshtools.tetgenstructs import NodePoint

from tetmeshtools.grid import Grid

def make_progress_test(end_x, end_y, end_z, steps=10, start_x=0, start_y=0, start_z=0):
    """
    make a progress test function
    Args:
        end_x (int): number of x columns
        end_y (int): number of y columns
        end_z (int): number of z columns
        steps (int): the number of iterations between reporting
        start_x (int): the count (zero start) of the first x column
        start_y (int): the count (zero start) of the first y column
        start_z (int): the count (zero start) of the first z column
    Returns:
        function(int, int, int) =>int or None
        (int): the total number of iterations
    """
    total = (end_x - start_x) * (end_y - start_y) * (end_z - start_z)
    stage = int(total/steps)

    def progress_test(current_iteration):
        """
        test if the current iteration should be reported for prograss
        Args:
            current_iteration (int): the number of the current iteration
        Returns:
            if the current iteration should be reported, the iteration , else
        """
        if current_iteration%stage == 0:
            print(f"Completed {current_iteration} out of {total} iterations")

    return progress_test

def make_density_at_vertex_grid(image, counts=None, progress=True):
    """
    Converts image into voxels with densities at vertices.
    Args:
        image (MRCImage): the input file
        counts ([int, int, int]): voxel counts on x, y and z axis, if not
        progress (bool): if true print out progress
    Returns:
        (Grid)
    """
    # get the start and end value of the image cube axis
    start = [image.x_origin,
             image.y_origin,
             image.z_origin]
    end = [image.cell_size[0]+start[0],
           image.cell_size[1]+start[1],
           image.cell_size[2]+start[2]]

    image_counts =[image.get_nx(), image.get_ny(), image.get_nz()]

    if counts is None:
        counts = [x-1 for x in image_counts]

    grid = Grid(counts, start, end, image_counts)
    if progress:
        print("Grid object constructed.", file=sys.stdout)
        print("Making vertices ...", file=sys.stdout)

    grid.build_grid(image)

    if progress:
        print("Vertices constructed", file=sys.stdout)

    return grid

def all_voxels_to_tets(image, counts, progress, use_six_tets=True):
    """
    Converts image into voxels of 5 tetrohedrons.
    Args:
        image (MRCImage): the input file
        counts ([int, int, int]): voxel counts on x, y and z axis
        progress (bool): if true print out progress
        five_tets (bool): it True use 5 tet decomp, else 6
    Returns:
        (Grid)
    """
    grid = make_density_at_vertex_grid(image, counts, progress)

    if use_six_tets:
        grid.build_six_tets()
    else:
        grid.build_five_tets()
    if progress:
        print("Tets constructed", file=sys.stdout)

    return grid

def convert_mrc_to_tets(input_file,
                        output_file,
                        threshold,
                        tetgen_out,
                        vtk_out,
                        verbose,
                        progress,
                        vox_counts,
                        low_vertices,
                        use_six_tets = False):
    """
    Converts the contents of a mrc file to a tetrohedron array (5 pre voxel).
    Args:
        input_file (pathlib.Path): name of input file
        output_file (pathlib.Path): name stem for output files
        threshold (float): the threshold below which results are ignored (isosurface value)
        tetgen_out (bool): if true produce tetgen format files
        vtk_out (bool): if true produce vtk file
        verbose (bool): if true give details of results
        progress (bool): if true print out progress
        vox_counts ([int]): voxels counts on x, y, z, if None use image voxel counts
        low_vertices (int): number of vertices below isovalue that triggers culling
        use_six_tets (bool): if True convert grid voxels to six tets, else five
    Returns:
        None
    """
    prune_level = v2t.PruneLevel(low_vertices)

    # start time
    time_start = process_time()

    with mrcfile.mmap(input_file, mode='r+') as mrc:
        image = mi.MRCImage(mrc)

        if vox_counts is None:
            vox_counts = [image.get_nx(), image.get_ny(), image.get_nz()]

        grid = all_voxels_to_tets(image, vox_counts, progress, use_six_tets)
        if progress:
            print(f"Grid formed {len(grid.get_vertices())} vertices")

        grid.crop_mesh_to_isovalue(threshold, prune_level, progress)
        if progress:
            print(f"Grid formed {len(grid.get_connectivities())} tets")

        grid.remove_surplas_vertices()
        if progress:
            print(f"redundant vertices removed {len(grid.get_vertices())} remain")

        if grid.get_total_num_voxels() <= 0:
            print(f"Error: threshold value of {threshold} yielded no results", file=sys.stderr)
            sys.exit()

        # end time
        time_end = process_time()

        if progress:
            print(f"Conversion in {time_end - time_start} seconds, writing files")

        if len(grid.get_connectivities()) == 0:
            print("No tetrahedrons were made, so no files written", file=sys.stdout)
            return

        if verbose:
            verbose_output(mrc,
                           grid,
                           use_six_tets)

        v2t.write_tets_to_files(grid.get_vertices(),
                                grid.get_connectivities(),
                                output_file,
                                tetgen_out,
                                vtk_out)

def verbose_output(mrc, grid, six_tets):
    """
    print description of the output
    Args:
        mrc (mrcfile): source file
        grid (Grid): mesh object
        six_tets (bool): if true gird uses 6 tetrahedra per voxel
    """
    print("\nSummary\n============")

    print(f"Number of voxels in mesh {grid.get_total_num_voxels()}")
    iv_size = voxel_size(mrc)
    message  = f"Voxel size in MRC image ({iv_size.dx:.3f}, {iv_size.dy:.3f}, {iv_size.dz:.3f})"
    message += f", volume {iv_size.dx * iv_size.dy * iv_size.dz:.3f}"
    print(message)

    mv_size = grid.get_voxel_size()
    vol = mv_size.dx * mv_size.dy * mv_size.dz
    print(f"Voxel size in mesh ({mv_size.dx:.3f}, {mv_size.dy:.3f}, {mv_size.dz:3f}), volume {vol:.3f}")

    if six_tets:
        print("6 tets per voxel")
        print(f"Number of tets with volume of {vol/6:.3f} is {len(grid.get_connectivities())}")
        return

    print("5 tets per voxel")
    verbose_output_five_tets(grid, vol)

def verbose_output_five_tets(grid, vol):
    """
    output the number of tets in each size
    Args:
        grid (Grid): mesh
        vol (float): volume of one mesh voxel
    """
    large = 0
    # to seperate tets vol/3 from vol/6 use vol/4
    limit = vol/4.0

    for tet in grid.get_connectivities():
        geom = []
        for index in tet:
            vert = grid.get_vertices()[index]
            geom.append(NodePoint(index, vert[0], vert[1], vert[2]))

        if tet_volume(geom) > limit:
            large += 1

    print(f"Number of tets with volume of {vol/3.0:.3f} is {large}")
    print(f"Number of tets with volume of {vol/6.0:.3f} is {len(grid.get_connectivities())-large}")
