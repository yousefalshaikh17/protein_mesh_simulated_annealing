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

import argparse
import pathlib
import mrcfile

from tetmeshtools.mrc_utility import voxel_size

def get_args():
    """
    get the command line arguments
        Returns
            (argparse.namespace)
    """
    parser = argparse.ArgumentParser(description="print an MRC file's voxel size")

    parser.add_argument("-i",
                        "--input",
                        type=pathlib.Path,
                        required=True,
                        help="input file")

    return parser.parse_args()

def main():
    """
    run the script
    """
    args = get_args()
    with mrcfile.open(args.input, mode='r+') as mrc:
        delta = voxel_size(mrc)

    print(f"{args.input} voxel size {round(delta.dx, 2)}, "\
          f"{round(delta.dy, 2)}, {round(delta.dz, 2)}")

if __name__ == "__main__":
    main()
