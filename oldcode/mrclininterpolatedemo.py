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
    tolerance = 0.00001
    origin = (-117.703995, 172.448, -120.83199)

    _, dist = image.density_or_distance_at(origin[0], origin[1]+10.0, origin[2]+10.0)
    if math.isclose(dist, 6.25, abs_tol=tolerance):
        print("PASS: x too low")
    else:
        print(f"FAIL: x too low ")

    _, dist = image.density_or_distance_at(origin[0]+500.0, origin[1]+10.0, origin[2]+10.0)
    if math.isclose(dist, 56.24999, abs_tol=tolerance):
        print("PASS: x too high")
    else:
        print(f"FAIL: x too high ")

    _, dist = image.density_or_distance_at(origin[0]+10.0, origin[1]+0.0, origin[2]+10.0)
    if math.isclose(dist, 6.25, abs_tol=tolerance):
        print("PASS: y too low")
    else:
        print("FAIL: y too low")


    _, dist = image.density_or_distance_at(origin[0]+10.0, origin[1]+500.0, origin[2]+10.0)
    if math.isclose(dist, 31506.250693, abs_tol=tolerance):
        print("PASS: y too high")
    else:
        print("FAIL: y too high unknown exception")


    _, dist = image.density_or_distance_at(origin[0]+10.0, origin[1]+10.0, origin[2]+0.0)
    if math.isclose(dist, 6.25, abs_tol=tolerance):
        print("PASS: z too low")
    else:
        print("FAIL: z too low ")

    _, dist = image.density_or_distance_at(origin[0]+10.0, origin[1]+10.0, origin[2]+500.0)
    if math.isclose(dist, 20306.25088, abs_tol=tolerance):
        print("PASS: z too high")
    else:
        print("FAIL: z too high")

def test_x(image, tolerance):
    """
    test interpolation on x axis
    Args:
        image (DensityTable)
        tolerance (float): epsilon for test
    """
    origin = (-117.703995, 172.448, -120.83199)
    value, _ = image.density_or_distance_at(origin[0]+100.0, origin[1]+100.0, origin[2]+100.0)
    target = image.data_at(20, 20, 20)
    if math.isclose(value, target, abs_tol=tolerance):
        print(f"TEST x pass (target {target})")
    else:
        print(f"TEST x fail: was {value} should be {target})")

    value, _ = image.density_or_distance_at(origin[0]+105.0, origin[1]+100.0, origin[2]+100.0)
    target = image.data_at(21, 20, 20)
    if math.isclose(value, target, abs_tol=tolerance):
        print(f"TEST x pass (target {target}")
    else:
        print(f"TEST x fail: was {value} should be {target}")

    value, _ = image.density_or_distance_at(origin[0]+102.5, origin[1]+100.0, origin[2]+100.0)
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
    origin = (-117.703995, 172.448, -120.83199)
    value, _ = image.density_or_distance_at(origin[0]+100.0, origin[1]+100.0, origin[2]+100.0)
    target = image.data_at(20, 20, 20)
    if math.isclose(value, target, abs_tol=tolerance):
        print(f"TEST y pass (target {target})")
    else:
        print(f"TEST y fail: was {value} should be {target})")

    value, _ = image.density_or_distance_at(origin[0]+100.0, origin[1]+105.0, origin[2]+100.0)
    target = image.data_at(20, 21, 20)
    if math.isclose(value, target, abs_tol=tolerance):
        print(f"TEST y pass (target {target}")
    else:
        print(f"TEST y fail: was {value} should be {target}")

    value, _ = image.density_or_distance_at(origin[0]+100.0, origin[1]+102.5, origin[2]+100.0)
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
    origin = (-117.703995, 172.448, -120.83199)
    value, _ = image.density_or_distance_at(origin[0]+100.0, origin[1]+100.0, origin[2]+100.0)
    target = image.data_at(20, 20, 20)
    if math.isclose(value, target, abs_tol=tolerance):
        print(f"TEST z pass (target {target})")
    else:
        print(f"TEST z fail: was {value} should be {target})")

    value, _ = image.density_or_distance_at(origin[0]+100.0, origin[1]+100.0, origin[2]+105.0)
    target = image.data_at(20, 20, 21)
    if math.isclose(value, target, abs_tol=tolerance):
        print(f"TEST z pass (target {target}")
    else:
        print(f"TEST z fail: was {value} should be {target}")

    value, _ = image.density_or_distance_at(origin[0]+100.0, origin[1]+100.0, origin[2]+102.5)
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
    tolerance = 0.000001

    if test == TestType.ALL:
        test_range(image)
        test_x(image, tolerance)
        test_y(image, tolerance)
        test_z(image, tolerance)
    elif test == TestType.RANGE:
        test_range(image)
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
    with mrcfile.open(args.input, mode='r+') as mrc:
        test_mrc(mrc, args.test)

if __name__ == "__main__":
    main()
