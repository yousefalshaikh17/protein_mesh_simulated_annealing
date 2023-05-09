import argparse
import pathlib
import numpy as np
import mrcfile
import itertools
from PIL import Image

import ffeamesh.mrclininterpolate as mli

def display_cross_section(image, xsize, zsize, y, isovalue, file_name):
    """
    display a section of the file
    Args:
        mrc (mrcfile): image to analyse
        args (argparse.namespace): options
    """
    x_range = np.linspace(image.offset_x, image.inner_size_x, xsize)
    z_range = np.linspace(image.offset_z, image.inner_size_z, zsize)

    slice = np.zeros((len(x_range), len(z_range)))

    for i, x in enumerate(x_range):
        for j, z in enumerate(z_range):
            slice[i, j] = image.density_at(x, y, z)

    show_image(slice, isovalue, file_name)

def show_image(slice_array, isovalue, file_name):
    image_array = np.zeros(slice_array.shape, dtype=np.uint8)

    for i in range(image_array.shape[0]):
        for j in range(image_array.shape[1]):
            if slice_array[i, j] > isovalue:
                image_array[i, j] = 255

    img = Image.fromarray(image_array, "L")
    # Display the Numpy array as Image
    #img.show()
    img.save(file_name)

def get_args():
    """
    get command line
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("-i",
                        "--input",
                        type=pathlib.Path,
                        required=True,
                        help="MRC source file")

    parser.add_argument("-x",
                        "--xsize",
                        type=int,
                        required=True,
                        help="x steps")

    parser.add_argument("-z",
                        "--zsize",
                        type=int,
                        required=True,
                        help="z steps")

    parser.add_argument("-v",
                        "--isovalue",
                        type=float,
                        required=True,
                        help="value defining surface")

    return parser.parse_args()

def main():
    """
    demo of linear interpolation in data from mrc file:
    map_equilibrium_mystructrigor_15A_0p00202.mrc
    """
    args = get_args()
    #with mrcfile.open(args.input, mode='r+') as mrc:
    with mrcfile.mrcmemmap.MrcMemmap(args.input) as mrc:
        image = mli.MRCImage(mrc)
        count = itertools.count(1)
        total = 305-45
        for i in range(45, 305, 1):
            y = float(i)
            file_name = f"yscan_{i}.png"
            display_cross_section(image, args.xsize, args.zsize, y, args.isovalue, file_name)
            print(f"Image {next(count)} of {total}")

if __name__ == "__main__":
    main()
