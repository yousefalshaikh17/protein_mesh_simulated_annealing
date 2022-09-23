# -*- coding: utf-8 -*-
#
#  This file is part of the FFEA simulation package
#
#  Copyright (c) by the Theory and Development FFEA teams,
#  as they appear in the README.md file.
#
#  FFEA is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  FFEA is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
#
#  To help us fund FFEA development, we humbly ask that you cite
#  the research papers on the package.
#

"""
Authors: Joanna Leng, Jonathan Pickering - University of Leeds
Emails: J.Leng@leeds.ac.uk, J.H.Pickering@leeds.ac.uk
"""
# set up linting conditions
# pylint: disable = import-error

from collections import namedtuple
import mrcfile
import numpy as np

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
    extract cell sizes from mrc file
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

def write_mrcfile(data, cell_props, out_file, label=None, overwrite=False):
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
            mrc.header.label[0] = label
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
        image

    Returns:
        (np.array): image with all voxels failing test set to zero
    """
    tmp = np.array([0 if value<threshold else value for value in image.flatten()], dtype=np.float32)
    return tmp.reshape(image.shape)
