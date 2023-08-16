import argparse
import pathlib
import mrcfile

import ffeamesh.mrclininterpolate as mi

def mrc_file(name):
    """
    check string is mrc image file
    Args:
        name (str)
    Returns
        pathlib.Path
    """
    path = pathlib.Path(name)

    if not path.suffix == ".mrc":
        raise ValueError()

    if not path.exists():
        raise ValueError()

    return path

def get_args():
    """
    get the command line arguments
        Returns
            (argparse.namespace)
    """
    parser = argparse.ArgumentParser("""test isosurface""")

    parser.add_argument("-i",
                        "--image_file",
                        type=mrc_file,
                        required=True,
                        help="mrc image file")

    return parser.parse_args()

def main():
    args = get_args()
    with mrcfile.open(args.image_file, mode='r+') as mrc:
        image = mi.MRCImage(mrc)
        for coord in range(0, 27):
            x_coord = float(coord)
            found_density, dist2 = image.density_or_distance_at(x_coord, 12.5, 12.5)
            print(f"{x_coord}, {found_density}, {dist2}")

if __name__ == "__main__":
    main()
