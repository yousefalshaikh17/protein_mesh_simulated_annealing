"""
storage, access and writer for the data in a Medical Research Council (MRC) format
file (https://www.ccpem.ac.uk/mrc_format/mrc2014.php). Reading is provided by
package 'mrcfile'.

-------------------------------------

Licensed under the GNU General Public License, Version 3.0 (the "License"); you
may not use this file except in compliance with the License. You may obtain a
copy of the License at <https://www.gnu.org/licenses/gpl-3.0.html>.

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

This work was funded by Joanna Leng's EPSRC funded RSE Fellowship (EP/R025819/1)

@copyright 2020
@author: j.h.pickering@leeds.ac.uk and j.leng@leeds.ac.uk
"""
# set up linting conditions
# pylint: disable = import-error

from collections import namedtuple
import mrcfile
import numpy as np

from tetmeshtools.voxelsize import VoxelSize

## a data struct for the mrcfile cell size
CellSize = namedtuple("CellSize", "x, y, z")

## a data struct for the mrcfile cell angles (degrees)
CellAngles = namedtuple("CellAngles", 'alpha, beta, gamma')

## a data struct for cell props
CellProps = namedtuple("CellProps", "size, angles")

def get_cell_props(mrc):
    """
    extract call size and angles to a CellProps object
    Args:
        mrc (mrcfile) input file
    Returns
        (CellProps) the cell size and angles of the input file
    """
    return CellProps(get_cell_sizes(mrc), get_cell_angles(mrc))

def get_cell_sizes(mrc):
    """
    extract cell sizes from mrc file  (voxel size x/nx etc)
    Args:
        mrc (mrcfile): the input file
    Returns
        (CellSize) the cell sizes from the input file
    """
    return CellSize(mrc.header.cella.x, mrc.header.cella.y, mrc.header.cella.z)

def get_cell_angles(mrc):
    """
    extract cell angles from mrc file
    Args:
        mrc (mrcfile): the input file
    Returns
        (CellAngles) the cell angles from the input file
    """
    return CellAngles(mrc.header.cellb.alpha,
                      mrc.header.cellb.beta,
                      mrc.header.cellb.gamma)

def write_mrcfile(data, cell_props, out_file, label=None, overwrite=False, origin=None):
    """
    write a numpy array to an mrc_file
    Args:
        data (numpy.ndarray) the image
        cell_size (CellSize): the unit cell dimensions
        out_file (pathlib.Path) the output file
        label (str): the header label for the file
        overwrite (bool): if true overwrite
    """
    if out_file.exists():
        if overwrite:
            out_file.unlink()
        else:
            print(f"File {out_file} exists, cannot overwrite!")
            return

    with mrcfile.mmap(out_file, mode='w+') as mrc:
        mrc.set_data(data)
        mrc.update_header_from_data()

        if label is not None:
            if isinstance(label, str):
                mrc.header.label[0] = label
            elif isinstance(label, list):
                for index in range(min(8, len(label))):
                    mrc.header.label[index] = label[index]

        if origin is not None:
            mrc.header.origin.x = origin[0]
            mrc.header.origin.y = origin[1]
            mrc.header.origin.z = origin[2]

        mrc.header.cella.x = cell_props.size.x
        mrc.header.cella.y = cell_props.size.y
        mrc.header.cella.z = cell_props.size.z
        mrc.header.cellb.alpha = cell_props.angles.alpha
        mrc.header.cellb.beta = cell_props.angles.beta
        mrc.header.cellb.gamma = cell_props.angles.gamma
        mrc.flush()

def threshold_mrc_image(image, threshold):
    """
    make a copy of  3D image with all voxels that have a density less
    than the threshold set to zero

    Args:
        image (np.array): the input array
        threshold (float): the limit: if value < thershold then value = 0.0

    Returns:
        (np.array): image with all voxels failing test set to zero
    """
    tmp = np.array([0 if value<threshold else value for value in image.flatten()], dtype=np.float32)
    return tmp.reshape(image.shape)

def voxel_size(mrc):
    """
    cacluate the voxel size
    Args:
        mrc (mrcfile): source data
    Returns
        VoxelSize
    """
    delta_x = mrc.header.cella.x/mrc.header.mx
    delta_y = mrc.header.cella.y/mrc.header.my
    delta_z = mrc.header.cella.z/mrc.header.mz
    return VoxelSize(delta_x, delta_y, delta_z)
