#!/usr/bin/env python3
"""
 A script that reads a mrc file and prints out its image intensity stats.

 -----------------------

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
import sys
from scipy import stats
import numpy as np
import mrcfile

def get_args():
    """
    get the command line arguments
        Returns
            (argparse.namespace)
    """
    description = ("reads a mrc file and prints out its "
                   "image intensity stats")
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("-i",
                        "--input",
                        type=pathlib.Path,
                        required=True,
                        help="input file")

    parser.add_argument("-o",
                        "--output",
                        type=pathlib.Path,
                        required=False,
                        help="output file")

    parser.add_argument("-w",
                        "--overwrite",
                        action="store_true",
                        help="overwrite")

    return parser.parse_args()

def validate_args(args):
    """
    check the command line arguments are sane
    Args
        args (argparse.Namespace)
    Returns:
        (str/None) error message if problem else None
    """
    if not args.input.exists():
        return f"input file {args.input} does not exist!"

    if args.output is not None and args.output.exists():
        if not args.overwrite:
            return f"output file {args.output} already exists run with -w to overwrite"

    return None

def get_image_stats(mrc):
    """
    generate the stats of the image array
    Args:
        mrc (mrcfile): the image source
    Returns:
        (scipy.stats.DescribeResult)
        (numpy.histogram)
    """
    descriptive_stats = stats.describe(mrc.data.flatten())
    histogram = np.histogram(mrc.data.flatten(), bins=10)

    return descriptive_stats, histogram

def print_image_stats(descriptive_stats,  histo, infile, outfile=None):
    """
    generate and print the stats of the image array
    Args:
        descriptive_stats (scipy.stats.DescribeResult): mean etc
        histo (numpy.histogram): ten bin histogram
        mrc (mrcfile): the image source
        infile (str): file name of input
        outfile (pathlib.Path) the output file
    """

    print(f"File, {infile}", file=outfile)
    print(f"Number of voxels in image, {descriptive_stats.nobs}", file=outfile)
    print(f"Minimum, {descriptive_stats.minmax[0]:.6}", file=outfile)
    print(f"Maximum, {descriptive_stats.minmax[1]:.6}", file=outfile)
    print(f"Mean, {descriptive_stats.mean:.6}", file=outfile)
    print(f"Variance {descriptive_stats.variance:.6}", file=outfile)

    print("\nbin, range, voxel count, fraction", file=outfile)
    for i, count in enumerate(histo[0]):
        fraction = count/descriptive_stats.nobs
        print(f"{i}, ({histo[1][i]:.6}, {histo[1][i+1]:.6}), {count}, {fraction:.3}", file=outfile)

def main():
    """
    run the script
    """
    args = get_args()
    message = validate_args(args)
    if message is not None:
        print(f"Error: {message}", file=sys.stderr, flush=True)
        sys.exit()

    with mrcfile.mmap(args.input, 'r') as mrc:
        descriptive_stats,  histogram = get_image_stats(mrc)
        if args.output is not None:
            with args.output.open("w") as outfile:
                print_image_stats(descriptive_stats,  histogram, str(args.input), outfile)
        else:
            print_image_stats(descriptive_stats,  histogram, str(args.input))

if __name__ == "__main__":
    main()
