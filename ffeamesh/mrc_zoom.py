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
        mrc_zoom.py
        Authors: Molly Gravett, Joanna Leng, Jarvellis Rogers - University of Leeds
        Emails: bsmgr@leeds.ac.uk, J.Leng@leeds.ac.uk, J.F.Rogers1@leeds.ac.uk
"""

# Zooms and rescales MRC files for tet_from_pix

from os import path
from scipy import ndimage
import mrcfile
import numpy as np
import scipy
import getopt
import sys

def mrc_zoom(mrcin, mrcdata, nx, ny, nz, resolution, mrcout):
    scale_factor = (mrcin.voxel_size.tolist()[0]) / resolution
    b = ndimage.zoom(mrcdata, scale_factor, order=3)
    mrcfilename_new = mrcout + ".mrc"

    with mrcfile.new(mrcfilename_new, overwrite=True) as mrc_2:
        mrc_2.set_data(b)
        mrc_2.header.origin = mrcin.header.origin
        old_vox = mrcin.voxel_size.tolist()
        new_list =[]
        nx_2 = mrc_2.header.nx
        ny_2 = mrc_2.header.ny
        nz_2 = mrc_2.header.nz
        scale_x = nx_2/nx
        scale_y = ny_2/ny
        scale_z = nz_2/nz
        new_list.append((1/scale_x)*old_vox[0]*(nx_2/(nx_2-1)))
        new_list.append((1/scale_y)*old_vox[1]*(ny_2/(ny_2-1)))
        new_list.append((1/scale_z)*old_vox[2]*(nz_2/(nz_2-1)))
        mrc_2.voxel_size = tuple(new_list)

    mrcfile.validate(mrcfilename_new)
    return mrcfilename_new

if __name__ == "__main__":
    def pathCheck(fpath):
        if not path.exists(fpath):
            msg = "ERROR: Input file path does not exist. Please try again."
            sys.exit(msg)

    def usage():
        helpMessage = """   Coding:   Molly Gravett (bsmgr@leeds.ac.uk), Joanna Leng (J.Leng@leeds.ac.uk), Jarvellis Rogers (J.F.Rogers1@leeds.ac.uk)

    mrc_zoom coarsens MRC files to a user-defined resolution and outputs them as a new .mrc file.

    General Options:
      -h [ --help ]             Print usage message.
      -i [ --input ] arg        File path to input MRC file. (required)
      -o [ --output ] arg       Name or /path/to/name of output MRC file. (required)
      -r [ --resolution ] arg   Resolution to coarsen input MRC file by. (required)
    """

        print(helpMessage)
        sys.exit()

    try:
        options, remainder = getopt.getopt(sys.argv[1:], "i:o:r:h", ["input=", "output=", "resolution=", "help"])
    except getopt.GetoptError as err:
        print("ERROR: " + str(err) + "\n")
        usage()

    for opt, arg in options:
        if opt in ("-i", "--input"):
            mrcfilename = arg
            pathCheck(mrcfilename)
        elif opt in ("-o", "--output"):
            chosen_filename = arg
        elif opt in ("-r", "--resolution"):
            try:
                resolution = float(arg)
            except ValueError:
                msg = "ERROR: Resolution must be a number."
                sys.exit(msg)
        elif opt in ("-h", "--help"):
            usage()

    # Checks to see if mandatory options have been called
    try:
        mrcfilename
    except NameError:
        msg = "ERROR: Input MRC file path not defined."
        sys.exit(msg)
    try:
        chosen_filename
    except NameError:
        msg = "ERROR: Output destination file path not defined."
        sys.exit(msg)
    try:
        resolution
    except NameError:
        msg = "ERROR: Resolution not defined."
        sys.exit(msg)

    mrc = mrcfile.open(mrcfilename, mode='r+')
    a = mrc.data
    nx = mrc.header.nx
    ny = mrc.header.ny
    nz = mrc.header.nz

    mrc_zoom(mrc, a, nx, ny, nz, resolution, chosen_filename)
