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
import pathlib
from ffeamesh.mrc_utility import (CellSize, CellAngles, CellProps, write_mrcfile)

def fill_test(image):
    """
    fill the image array
    Args:
        image: numpy.ndarray
    """
    for i in range(1, 4):
        for j in range(1, 4):
            for k in range(1, 4):
                image[i][j][k] = 0.8

    image[2][2][2] = 1.2

def main():
    """
    script produces a 5x5x5 mrc file with the value 1.2 at the centre
    the value 0.8 around that and value 0.2 at the edge. The image is
    saved as file simple_test.mrc.
    """
    cell_size = 25.0
    cell_angle = 90.0
    label="Simple test"
    test_image = np.full((5, 5, 5), dtype=np.float32, fill_value=0.2)
    fill_test(test_image)

    props = CellProps(CellSize(cell_size, cell_size, cell_size),
                      CellAngles(cell_angle, cell_angle, cell_angle))

    write_mrcfile(test_image, props, pathlib.Path("simple_test.mrc"), label, True)

if __name__ == "__main__":
    main()
