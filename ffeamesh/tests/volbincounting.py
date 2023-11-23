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

import unittest
from ffeamesh.utility import VolumeBins

class TestVolBins(unittest.TestCase):
    """
    tests of Tetgen file reader
    """

    def setUp(self):
        """
        build a full test class
        """
        pass

    def tearDown(self):
        """
        clean up
        """
        pass

    def test_nodepoint_to_edge_array(self):
        """
        test NodePoint to_edge_arry
        """
        vol_bins = VolumeBins()

        for i in [1.0, 1.23, 1.7, 2.3, 1.1, 2.7]:
            vol_bins.add(i)

        self.assertEqual()

def main():
    """
    run a test
    """
    vol_bins = VolumeBins()

    for i in [1.0, 1.23, 1.7, 2.3, 1.1, 2.7]:
        vol_bins.add(i)

    print("Should be:")
    print("1.0 => 3")
    print("2.0 => 2")
    print("3.0 => 1")

    print("Was:")
    for bin, count in vol_bins.bin_counts.items():
        print(f"{bin} => {count}")

if __name__ == '__main__':
    main()
