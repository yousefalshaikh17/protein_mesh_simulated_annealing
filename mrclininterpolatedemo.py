import enum
import argparse
import pathlib
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

def test_x(image):
    """
    test interpolation on x axis
    Args:
        image (DensityTable)
    """
    print("test x")

def test_y(image):
    """
    test interpolation on y axis
    Args:
        image (DensityTable)
    """
    print("test y")

def test_z(image):
    """
    test interpolation on z axis
    Args:
        image (DensityTable)
    """
    print("test z")

def test_mrc(mrc, test):
    """
    select test to run
    Args:
        mrc (mrcfile): image to analyse
        test (TestType): test option
    """
    image = mli.MRCImage(mrc)

    if test == TestType.ALL:
        test_range(image)
        test_x(image)
        test_y(image)
        test_z(image)
        return
    elif test == TestType.RANGE:
        test_range(image)
        return
    elif test == TestType.INTERPX:
        test_x(image)
        return
    elif test == TestType.INTERPY:
        test_y(image)
        return
    elif test == TestType.INTERPZ:
        test_z(image)
        return
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
