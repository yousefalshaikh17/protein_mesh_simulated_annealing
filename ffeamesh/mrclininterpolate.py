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
import numpy as np
import collections
import mrcfile

import ffeamesh.mrc_utility as mu

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

        self.dx = mrc.header.cella.x/mrc.header.mx
        self.dy = mrc.header.cella.y/mrc.header.my
        self.dz = mrc.header.cella.z/mrc.header.mz

        self.offset_x = self.dx/2.0
        self.offset_y = self.dy/2.0
        self.offset_z = self.dz/2.0

        self.inner_size_x = self.dx*(self._file.header.nx-1)
        self.inner_size_y = self.dy*(self._file.header.ny-1)
        self.inner_size_z = self.dz*(self._file.header.nz-1)

    def test_inner_coords(self, x, y, z):
        """"
        ensure offset coordinates are in inner array range
        Args:
            x (float): inner x coordinate
            y (float): inner y coordinate
            z (float): inner z coordinate
        Raise
            ValueError if out of range
        """
        message = "coordinate below intercept calculation range"
        if x < 0.0:
            raise ValueError(f"x {message}")
        if y < 0.0:
            raise ValueError(f"y {message}")
        if z < 0.0:
            raise ValueError(f"z {message}")

        message = "coordinate above intercept calculation range"
        if x >= self.inner_size_x:
            raise ValueError(f"x {message}")
        if y >= self.inner_size_y:
            raise ValueError(f"y {message}")
        if z >= self.inner_size_z:
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

        x_frac, x_index = np.modf(offset_x/self.dx)
        y_frac, y_index = np.modf(offset_y/self.dy)
        z_frac, z_index = np.modf(offset_z/self.dz)

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

# class DummyMRC():
#     """
#     dummy for mrc file holder
#     """

#     def __init__(self):
#         self.nx = 5
#         self.ny = 5
#         self.nz = 5
#         self.dx = 0.3
#         self.dy = 0.4
#         self.dz = 0.5

#         self.offset_x = self.dx/2.0
#         self.offset_y = self.dy/2.0
#         self.offset_z = self.dz/2.0

#         self.inner_size_x = self.dx*(self.nx-1)
#         self.inner_size_y = self.dy*(self.ny-1)
#         self.inner_size_z = self.dz*(self.nz-1)

#         self.data = np.zeros((5, 5, 5))
#         for i in range(5):
#             for j in range(5):
#                 for k in range(5):
#                     self.data[i][j][k] = float(k+2*j+3*i)

    # def test_inner_coords(self, x, y, z):
    #     """"
    #     ensure offset coordinates are in inner array range
    #     Args:
    #         x (float): inner x coordinate
    #         y (float): inner y coordinate
    #         z (float): inner z coordinate
    #     Raise
    #         ValueError if out of range
    #     """
    #     if x < 0.0:
    #         raise ValueError("x coord under range")
    #     if y < 0.0:
    #         raise ValueError("y coord under range")
    #     if z < 0.0:
    #         raise ValueError("z coord under range")

    #     if x >= self.inner_size_x:
    #         raise ValueError("x coord over range")
    #     if y >= self.inner_size_y:
    #         raise ValueError("y coord over range")
    #     if z >= self.inner_size_z:
    #         raise ValueError("z coord over range")

    # def to_coords(self, image_x, image_y, image_z):
    #     """
    #     convert image coordinates to inner coordinates
    #     Args:
    #         image_x (float): image x coordinate
    #         image_y (float): image y coordinate
    #         image_z (float): image z coordinate
    #     Returns
    #         InnerCoords: the array index and frac
    #     """
    #     offset_x = image_x - self.offset_x
    #     offset_y = image_y - self.offset_y
    #     offset_z = image_z - self.offset_z

    #     self.test_inner_coords(offset_x, offset_y, offset_z)

    #     x_frac, x_index = np.modf(offset_x/self.dx)
    #     y_frac, y_index = np.modf(offset_y/self.dy)
    #     z_frac, z_index = np.modf(offset_z/self.dz)

    #     return make_inner_coords(round(x_index), x_frac,
    #                              round(y_index), y_frac,
    #                              round(z_index), z_frac,)

    # def linear_interp(self, coords):
    #     """
    #     find linear intrpolated density
    #     Args:
    #         coords (InnerCoords): coordinates
    #     Returns:
    #         (float): interpolated image denisty
    #     """
    #     a000 = self.data[coords.z_index, coords.y_index, coords.x_index]
    #     a100 = self.data[coords.z_index, coords.y_index, coords.x_index+1]
    #     a010 = self.data[coords.z_index, coords.y_index+1, coords.x_index]
    #     a001 = self.data[coords.z_index+1, coords.y_index, coords.x_index]
    #     a110 = self.data[coords.z_index, coords.y_index+1, coords.x_index+1]
    #     a101 = self.data[coords.z_index+1, coords.y_index, coords.x_index+1]
    #     a011 = self.data[coords.z_index+1, coords.y_index+1, coords.x_index]
    #     a111 = self.data[coords.z_index+1, coords.y_index+1, coords.x_index+1]

    #     value =  a000*coords.x_frac_comp*coords.y_frac_comp*coords.z_frac_comp
    #     value += a001*coords.x_frac_comp*coords.y_frac_comp*coords.z_frac
    #     value += a010*coords.x_frac_comp*coords.y_frac*coords.z_frac_comp
    #     value += a011*coords.x_frac_comp*coords.y_frac*coords.z_frac
    #     value += a100*coords.x_frac*coords.y_frac_comp*coords.z_frac_comp
    #     value += a101*coords.x_frac*coords.y_frac_comp*coords.z_frac
    #     value += a110*coords.x_frac*coords.y_frac*coords.z_frac_comp
    #     value += a111*coords.x_frac*coords.y_frac*coords.z_frac

    #     return value

    # def density_at(self, image_x, image_y, image_z):
    #     """
    #     return the interpolated density at image point
    #     Args:
    #         image_x (float): x image coordinate
    #         image_y (float): y image coordinate
    #         image_z (float): z image coordinate
    #     Retruns
    #         float: interpolated density
    #     Raises:
    #         ValueError if x, y or z out of range
    #     """
    #     coords = self.to_coords(image_x, image_y, image_z)
    #     return self.linear_interp(coords)