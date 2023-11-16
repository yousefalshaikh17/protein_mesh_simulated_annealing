"""
 tetvolume,py
 
 ----------------------------

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

from ffeamesh.fivetets import tet_volume

def main():
    """
    test using (0,0,0) (18,0,0) (9,18,0) (9,9,9)
    """
    target = 486.0

    coords = []
    coords.append([0.0, 0.0, 0.0])
    coords.append([18.0, 0.0, 0.0])
    coords.append([9.0, 18.0, 0.0])
    coords.append([9.0, 9.0, 9.0])

    volume = tet_volume(coords)
    message = "{}: volume calculated was {} should be {}."

    if abs(target-volume) > 0.0:
        print(message.format("Fail", volume, target))
    else:
        print(message.format("Pass", volume, target))

if __name__ == '__main__':
    main()
