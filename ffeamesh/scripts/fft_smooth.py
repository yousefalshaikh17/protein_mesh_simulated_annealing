#!/usr/bin/env python3
"""
 This file is part of the FFEA simulation package

 Copyright (c) by the Theory and Development FFEA teams,
 as they appear in the README.md file.

 FFEA is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 FFEA is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with FFEA.  If not, see <http://www.gnu.org/licenses/>.

 To help us fund FFEA development, we humbly ask that you cite
 the research papers on the package.

Created on Mon Oct 25 11:50:39 2021

@author: mollygravett
modified jonathan pickering 23Aug22
"""
# set up linting conditions
# pylint: disable = import-error
import argparse
import pathlib
from enum import Enum
import mrcfile
from scipy import ndimage
import numpy as np

class Algorithm(Enum):
    """
    enumeration of possible kernels
    """
    ELLIPSOID = 'ellipsoid'
    GAUSSIAN = 'gaussian'
    UNIFORM = 'uniform'

def get_args():
    """
    get the command line arguments
        Returns
            (argparse.namespace)
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("-i",
                        "--input",
                        type=pathlib.Path,
                        required=True,
                        help="input file")

    parser.add_argument("-o",
                        "--output",
                        type=pathlib.Path,
                        required=True,
                        help="output file")

    parser.add_argument("-p",
                        "--parameter",
                        type=float,
                        default=5.0,
                        help="""size of sphere or box used in ellipsoid or uniform
                            filtering; or the std dev for gaussian filtering""")

    parser.add_argument("-a",
                        "--algorithm",
                        type=Algorithm,
                        default=Algorithm.ELLIPSOID,
                        help="choice of algorithm used in fft filtering",
                        choices=list(Algorithm))

    return parser.parse_args()

def refine_volume_data_fft(in_file, out_file, param, algorithm):
    """
    simplify and smooth volume data in file
        Args:
            in_file (pathlib.Path) the input file
            out_file (pathlib.Path) the output file
            param (float) size of the ellipsoid/box or sigma used in filtering
            algorithm (Algorithm): the type of fft to be used
    """
    with mrcfile.open(in_file, mode='r+') as mrc:
        tmp = np.fft.fftn(mrc.data)

        if algorithm == Algorithm.ELLIPSOID:
            tmp = ndimage.fourier_ellipsoid(tmp, size=param)
        elif algorithm == Algorithm.GAUSSIAN:
            tmp = ndimage.fourier_gaussian(tmp, sigma=param)
        elif algorithm == Algorithm.UNIFORM:
            tmp = ndimage.fourier_uniform(tmp, size=param)

        result = np.fft.ifftn(tmp)
        result = result.real
        result = result.astype(np.single)

        output_mrcfile(result, out_file, mrc)

def output_mrcfile(data, out_file, original_mrc):
    """
    write out a scaled mrc file
        Args:
            data (numpy.ndarray): the volume data
            out_file (pathlib.Path) the output file
            original_mrc (mrcfile.File) the original input file
    """
    with mrcfile.new(out_file, overwrite=True) as mrc:
        mrc.set_data(data)
        mrc.header.origin = original_mrc.header.origin
        mrc.voxel_size = original_mrc.voxel_size

def main():
    """
    run the script
    """
    args = get_args()
    refine_volume_data_fft(args.input, args.output, args.parameter, args.algorithm)
    print(f"{args.input} filtered and written to {args.output}")

if __name__ == "__main__":
    main()
