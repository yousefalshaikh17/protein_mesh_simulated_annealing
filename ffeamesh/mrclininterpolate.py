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

        self.delta_x = mrc.header.cella.x/mrc.header.mx
        self.delta_y = mrc.header.cella.y/mrc.header.my
        self.delta_z = mrc.header.cella.z/mrc.header.mz

        self.offset_x = self.delta_x/2.0
        self.offset_y = self.delta_y/2.0
        self.offset_z = self.delta_z/2.0

        self.inner_size_x = self.delta_x*(self._file.header.nx-1)
        self.inner_size_y = self.delta_y*(self._file.header.ny-1)
        self.inner_size_z = self.delta_z*(self._file.header.nz-1)

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
        message = "coordinate below intercept calculation range"
        if x_coord < 0.0:
            raise ValueError(f"x {message}")
        if y_coord < 0.0:
            raise ValueError(f"y {message}")
        if z_coord < 0.0:
            raise ValueError(f"z {message}")

        message = "coordinate above intercept calculation range"
        if x_coord >= self.inner_size_x:
            raise ValueError(f"x {message}")
        if y_coord >= self.inner_size_y:
            raise ValueError(f"y {message}")
        if z_coord >= self.inner_size_z:
            raise ValueError(f"z {message}")

    def to_coords(self, image_x, image_y, image_z):
        """
        convert image coordinates to inner coordinates
        Args:
            image_x (float): image x coordinate
            image_y (float): image y coordinate
            image_z (float): image z coordinate
        Returns
            InnerCoords: the array index and frac
        """
        offset_x = image_x - self.offset_x
        offset_y = image_y - self.offset_y
        offset_z = image_z - self.offset_z

        self.test_inner_coords(offset_x, offset_y, offset_z)

        x_frac, x_index = np.modf(offset_x/self.delta_x)
        y_frac, y_index = np.modf(offset_y/self.delta_y)
        z_frac, z_index = np.modf(offset_z/self.delta_z)

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

    def density_at(self, image_x, image_y, image_z):
        """
        return the interpolated density at image point
        Args:
            image_x (float): x image coordinate
            image_y (float): y image coordinate
            image_z (float): z image coordinate
        Retruns
            float: interpolated density
        Raises:
            ValueError if x, y or z out of range
        """
        coords = self.to_coords(image_x, image_y, image_z)
        return self.linear_interp(coords)

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
