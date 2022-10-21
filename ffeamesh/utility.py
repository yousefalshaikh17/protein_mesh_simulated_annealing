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

    Authors: Joanna Leng, Jonathan Pickering - University of Leeds
    Emails: J.Leng@leeds.ac.uk, J.H.Pickering@leeds.ac.uk
"""
# set up linting
# pylint: disable = import-error

import numpy as np


def tet_volume(coords):
    """
    Find tet volume by ((side1 x side2).side2)/6
    Args:
        coords [CoordTransform x 4]: the vertices of the tet
    Returns:
        (float): the volume of the tet
    """
    # inner function
    def coords_to_np_vec(start, end):
        """
        make vectors (arrays) from cartesian coordinates of two CoordTransform
        Args:
            start (float list): start point
            end (float list): end point
        Returns:
            [x, y, z]: vector from start to end
        """
        vx = end[0] - start[0]
        vy = end[1] - start[1]
        vz = end[2] - start[2]

        return [vx, vy, vz]

    sides = []
    for coord in coords[1:]:
        sides.append(coords_to_np_vec(coords[0], coord))

    return abs(np.dot(np.cross(sides[0], sides[1]), sides[2]))/6.0
    
    
    
class VolumeBins():
    """
    Store for histogram of tet volumes every integer is a bin
    TODO Significant digits.
    """
    def __init__(self):
        """
        initalize the object
        """
        ## the bins
        self.bin_counts = {}

    def __repr__(self):
        """
        return string rep of object
        """
        return f"<VolumeBins : bins {len(self.bin_counts)}>"

    def add(self, vol):
        """
        add a volume to the bin count
        Args:
            vol (float): the volume to be added
        """
        tmp = round(vol, 0)
        if tmp in self.bin_counts.keys():
            self.bin_counts[tmp] += 1
        else:
            self.bin_counts[tmp] = 1
            
            
def print_voxel_stats(points, tet_connectities):
    """
    Print stats on volumes (as a text histogram).
    """
    vol_histo = VolumeBins()

    for tet in tet_connectities:
        coords = []
        for index in tet:
            coords.append(points[index])
        vol_histo.add(tet_volume(coords))

    for vol, count in vol_histo.bin_counts.items():
        print(f"Volume {vol}: {count} tets")
