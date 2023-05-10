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

import enum
import argparse
import pathlib
import math
import mrcfile

import ffeamesh.mrclininterpolate as mli

class TestType(enum.Enum):
    """
    possible tests
    """
    ALL     = "all"
    RANGE   = "range"
    INTERPX = "interpx"
    INTERPY = "interpy"
    INTERPZ = "interpz"

def test_range(image):
    """
    test allowed ranges
    Args:
        image (DummyMRC): the source
    """
    try:
        image.density_at(0.0, 10.0, 10.0)
        print("FAIL: x too low no exception")
    except ValueError as err:
        if str(err) == "x coord under range":
            print("PASS: x too low")
        else:
            print("FAIL: x too low unknown exception")

    try:
        image.density_at(500.0, 10.0, 10.0)
        print("FAIL: x too high no exception")
    except ValueError as err:
        if str(err) == "x coord over range":
            print("PASS: x too high")
        else:
            print(f"FAIL: x too high unknown exception {err}")

    try:
        image.density_at(10.0, 0.0, 10.0)
        print("FAIL: y too low no exception")
    except ValueError as err:
        if str(err) == "y coord under range":
            print("PASS: y too low")
        else:
            print("FAIL: y too low unknown exception")

    try:
        image.density_at(10.0, 500.0, 10.0)
        print("FAIL: y too high no exception")
    except ValueError as err:
        if str(err) == "y coord over range":
            print("PASS: y too high")
        else:
            print("FAIL: y too high unknown exception")

    try:
        image.density_at(10.0, 10.0, 0.0)
        print("FAIL: z too low no exception")
    except ValueError as err:
        if str(err) == "z coord under range":
            print("PASS: z too low")
        else:
            print("FAIL: z too low unknown exception")

    try:
        image.density_at(10.0, 10.0, 500.0)
        print("FAIL: z too high no exception")
    except ValueError as err:
        if str(err) == "z coord over range":
            print("PASS: z too high")
        else:
            print("FAIL: z too high unknown exception")

def test_x(image, tolerance):
    """
    test interpolation on x axis
    Args:
        image (DensityTable)
        tolerance (float): epsilon for test
    """
    value = image.density_at(102.5, 102.5, 102.5)
    target = image.data_at(20, 20, 20)
    if math.isclose(value, target, abs_tol=tolerance):
        print(f"TEST x pass (target {target})")
    else:
        print(f"TEST x fail: was {value} should be {target})")

    value = image.density_at(107.5, 102.5, 102.5)
    target = image.data_at(21, 20, 20)
    if math.isclose(value, target, abs_tol=tolerance):
        print(f"TEST x pass (target {target}")
    else:
        print(f"TEST x fail: was {value} should be {target}")

    value = image.density_at(105.0, 102.5, 102.5)
    target = image.data_at(21, 20, 20)*0.5 + image.data_at(20, 20, 20)*0.5
    if math.isclose(value, target, abs_tol=tolerance):
        print(f"TEST x pass (target {target})")
    else:
        print(f"TEST x fail: was {value} should be {target}")

def test_y(image, tolerance):
    """
    test interpolation on y axis
    Args:
        image (DensityTable)
        tolerance (float): epsilon for test
    """
    value = image.density_at(102.5, 102.5, 102.5)
    target = image.data_at(20, 20, 20)
    if math.isclose(value, target, abs_tol=tolerance):
        print(f"TEST y pass (target {target})")
    else:
        print(f"TEST y fail: was {value} should be {target})")

    value = image.density_at(102.5, 107.5, 102.5)
    target = image.data_at(20, 21, 20)
    if math.isclose(value, target, abs_tol=tolerance):
        print(f"TEST y pass (target {target}")
    else:
        print(f"TEST y fail: was {value} should be {target}")

    value = image.density_at(102.5, 105.0, 102.5)
    target = image.data_at(20, 21, 20)*0.5 + image.data_at(20, 20, 20)*0.5
    if math.isclose(value, target, abs_tol=tolerance):
        print(f"TEST y pass (target {target})")
    else:
        print(f"TEST y fail: was {value} should be {target}")

def test_z(image, tolerance):
    """
    test interpolation on z axis
    Args:
        image (DensityTable)
        tolerance (float): epsilon for test
    """
    value = image.density_at(102.5, 102.5, 102.5)
    target = image.data_at(20, 20, 20)
    if math.isclose(value, target, abs_tol=tolerance):
        print(f"TEST z pass (target {target})")
    else:
        print(f"TEST z fail: was {value} should be {target})")

    value = image.density_at(102.5, 102.5, 107.5)
    target = image.data_at(20, 20, 21)
    if math.isclose(value, target, abs_tol=tolerance):
        print(f"TEST z pass (target {target}")
    else:
        print(f"TEST z fail: was {value} should be {target}")

    value = image.density_at(102.5, 102.5, 105.0)
    target = image.data_at(20, 20, 21)*0.5 + image.data_at(20, 20, 20)*0.5
    if math.isclose(value, target, abs_tol=tolerance):
        print(f"TEST z pass (target {target})")
    else:
        print(f"TEST z fail: was {value} should be {target}")

def test_mrc(mrc, test):
    """
    select test to run
    Args:
        mrc (mrcfile): image to analyse
        test (TestType): test option
    """
    image = mli.MRCImage(mrc)
    tolerance = 0.00000000000001

    if test == TestType.ALL:
        test_range(image)
        test_x(image, tolerance)
        test_y(image, tolerance)
        test_z(image, tolerance)
    elif test == TestType.RANGE:
        test_range(image, tolerance)
    elif test == TestType.INTERPX:
        test_x(image, tolerance)
    elif test == TestType.INTERPY:
        test_y(image, tolerance)
    elif test == TestType.INTERPZ:
        test_z(image, tolerance)
    else:
        print(f"Error unknown test {test}")

def get_args():
    """
    get command line
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("-t",
                        "--test",
                        type=TestType,
                        required=True,
                        help = "test: " + str([e.value for e in TestType]))

    parser.add_argument("-i",
                        "--input",
                        type=pathlib.Path,
                        required=True,
                        help="MRC source file")

    return parser.parse_args()

def main():
    """
    run the demo
    """
    args = get_args()
    print(args)
    with mrcfile.open(args.input, mode='r+') as mrc:
        test_mrc(mrc, args.test)

if __name__ == "__main__":
    main()
