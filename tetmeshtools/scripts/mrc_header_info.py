#!/usr/bin/env python3
"""
 Script that makes simple mrc image files for use in testing.

 ------------------------

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

def get_args():
    """
    get the command line arguments
        Returns
            (argparse.namespace)
    """
    parser = argparse.ArgumentParser(description="print the header from an MRC file")

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

    if not args.input.exists():
        print(f"Error file {args.input} does not exist.")
        return

    with mrcfile.mmap(args.input, 'r') as mrc:
        for item in mrc.header.dtype.names:
            if item != 'label':
                print(f"{item} => {mrc.header[item]}")
            else:
                print(type(item))
                for index, label_s in enumerate(mrc.header[item]):
                    if label_s:
                        print(f'Label {index}:\t{label_s.decode()}')

if __name__ == "__main__":
    main()
