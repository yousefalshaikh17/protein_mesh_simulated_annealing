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
import collections
import numpy as np

## container for the index of cell and fractional coords inside cell
InnerCoords = collections.namedtuple("InnerCoords", ["x_index",
                                                     "x_frac",
                                                     "x_frac_comp",
                                                     "y_index",
                                                     "y_frac",
                                                     "y_frac_comp",
                                                     "z_index",
                                                     "z_frac",
                                                     "z_frac_comp"])

def make_inner_coords(x_index, x_frac, y_index, y_frac, z_index, z_frac):
    """
    make an inner coords object from indices and fractions
    Args:
        x_index (int): x array index
        x_frac (float): x fractional coordinate in voxel
        y_index (int): y array index
        y_frac (float): y fractional coordinate in voxel
        z_index (int): z array index
        z_frac (float): z fractional coordinate in voxel
    Retuns:
        InnerCoords
    """
    x_frac_comp = 1.0 - x_frac
    y_frac_comp = 1.0 - y_frac
    z_frac_comp = 1.0 - z_frac

    return InnerCoords(x_index, x_frac, x_frac_comp,
                       y_index, y_frac, y_frac_comp,
                       z_index, z_frac, z_frac_comp,)

class MRCImage():
    """
    container for mrcfile providing functions for
    linerar interpolation of image densities
    """

    def __init__(self, mrc):
        """
        set up object
        Args:#
            mrc (mrcfile.mrcfile): the  mrc file
        """

        ## the mrc file
        self._file = mrc

        # get coordinate origin
        self.x_origin = np.float32(mrc.header.origin['x'])
        self.y_origin = np.float32(mrc.header.origin['y'])
        self.z_origin = np.float32(mrc.header.origin['z'])

        self.delta_x = mrc.header.cella.x/mrc.header.mx
        self.delta_y = mrc.header.cella.y/mrc.header.my
        self.delta_z = mrc.header.cella.z/mrc.header.mz

        self.offset_x = self.delta_x/2.0
        self.offset_y = self.delta_y/2.0
        self.offset_z = self.delta_z/2.0

        self.inner_size_x = self.delta_x*(self._file.header.nx-1)
        self.inner_size_y = self.delta_y*(self._file.header.ny-1)
        self.inner_size_z = self.delta_z*(self._file.header.nz-1)

        self.low_limit_x = self.x_origin + self.offset_x
        self.low_limit_y = self.y_origin + self.offset_y
        self.low_limit_z = self.z_origin + self.offset_z

        self.high_limit_x = self.low_limit_x + self.inner_size_x
        self.high_limit_y = self.low_limit_y + self.inner_size_y
        self.high_limit_z = self.low_limit_z + self.inner_size_z

    def test_inner_coords(self, x_coord, y_coord, z_coord):
        """"
        ensure offset coordinates are in inner array range
        Args:
            x_coord (float): inner x coordinate
            y_coord (float): inner y coordinate
            z_coord (float): inner z coordinate
        Raise
            ValueError if out of range
        """
        # binary flag
        flag = 0
        if x_coord < self.low_limit_x:
            flag += 1
        if y_coord < self.low_limit_y:
            flag += 2
        if z_coord < self.low_limit_z:
            flag += 4

        if x_coord >= self.high_limit_x:
            flag += 8
        if y_coord >= self.high_limit_y:
            flag += 16
        if z_coord >= self.high_limit_z:
            flag += 32

        return flag

    def to_coords(self, image_x, image_y, image_z):
        """
        convert image coordinates to inner coordinates
        Args:
            image_x (float): image x coordinate
            image_y (float): image y coordinate
            image_z (float): image z coordinate
        Returns
            if (x, y, z) in image
                InnerCoords: the array index and frac
            else
                float distance from image cuboid
        """
        flag = self.test_inner_coords(image_x, image_y, image_z)
        if flag != 0:
            return flag

        true_x = image_x - self.x_origin
        true_y = image_y - self.y_origin
        true_z = image_z - self.z_origin

        x_frac, x_index = np.modf(true_x/self.delta_x)
        y_frac, y_index = np.modf(true_y/self.delta_y)
        z_frac, z_index = np.modf(true_z/self.delta_z)

        return make_inner_coords(round(x_index), x_frac,
                                 round(y_index), y_frac,
                                 round(z_index), z_frac,)

    def linear_interp(self, coords):
        """
        find linear intrpolated density
        Args:
            coords (InnerCoords): coordinates
        Returns:
            (float): interpolated image denisty
        """
        a000 = self._file.data[coords.z_index, coords.y_index, coords.x_index]
        a100 = self._file.data[coords.z_index, coords.y_index, coords.x_index+1]
        a010 = self._file.data[coords.z_index, coords.y_index+1, coords.x_index]
        a001 = self._file.data[coords.z_index+1, coords.y_index, coords.x_index]
        a110 = self._file.data[coords.z_index, coords.y_index+1, coords.x_index+1]
        a101 = self._file.data[coords.z_index+1, coords.y_index, coords.x_index+1]
        a011 = self._file.data[coords.z_index+1, coords.y_index+1, coords.x_index]
        a111 = self._file.data[coords.z_index+1, coords.y_index+1, coords.x_index+1]

        value =  a000*coords.x_frac_comp*coords.y_frac_comp*coords.z_frac_comp
        value += a001*coords.x_frac_comp*coords.y_frac_comp*coords.z_frac
        value += a010*coords.x_frac_comp*coords.y_frac*coords.z_frac_comp
        value += a011*coords.x_frac_comp*coords.y_frac*coords.z_frac
        value += a100*coords.x_frac*coords.y_frac_comp*coords.z_frac_comp
        value += a101*coords.x_frac*coords.y_frac_comp*coords.z_frac
        value += a110*coords.x_frac*coords.y_frac*coords.z_frac_comp
        value += a111*coords.x_frac*coords.y_frac*coords.z_frac

        return value

    def density_or_distance_at(self, image_x, image_y, image_z):
        """
        return the interpolated density at image point
        Args:
            image_x (float): x image coordinate
            image_y (float): y image coordinate
            image_z (float): z image coordinate
        Retruns
            (float, float): interpolated density and distance from image squared
        """
        coords = self.to_coords(image_x, image_y, image_z)

        if isinstance(coords, int):
            dist2 = self.dist_to_image_squared(coords, image_x, image_y, image_z)
            return 0.0, dist2

        density = self.linear_interp(coords)
        return density, 0.0

    def data_at(self, voxel_x, voxel_y, voxel_z):
        """
        get image density for voxel
        Args:
            voxel_x (int): array index
            voxel_y (int): array index
            voxel_z (int): array index
        Returns:
            float
        """
        return self._file.data[voxel_z, voxel_y, voxel_x]

    def dist_to_image_squared(self, region_flag, image_x, image_y, image_z):
        """
        return the square of the distance from the point to the image cuboid
        Args:
            region_flag (int): binary flag
            image_x (float): x coordinate in image coordinates
            image_y (float): y coordinate in image coordinates
            image_z (float): z coordinate in image coordinates
        returns:
            (float)
        """
        x_coord = 0.0
        y_coord = 0.0
        z_coord = 0.0

        # lower limits
        if region_flag & 1:
            x_coord = image_x - self.low_limit_x
        if region_flag & 2:
            y_coord = image_y - self.low_limit_y
        if region_flag & 4:
            z_coord = image_z - self.low_limit_z

        # high limits
        if region_flag & 8:
            x_coord = image_x - self.high_limit_x
        if region_flag & 16:
            y_coord = image_y - self.high_limit_y
        if region_flag & 32:
            z_coord = image_z - self.high_limit_z

        return x_coord*x_coord + y_coord*y_coord + z_coord*z_coord
