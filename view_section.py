import argparse
import pathlib
import numpy as np
import mrcfile
from PIL import Image

import ffeamesh.mrclininterpolate as mli

def display_cross_section(image, args):
    """
    display a section of the file
    Args:
        mrc (mrcfile): image to analyse
        args (argparse.namespace): options
    """
    x_range = np.linspace(image.offset_x, image.inner_size_x, args.xsize)
    y_range = np.linspace(image.offset_y, image.inner_size_y, args.xsize)
    z = args.zcoord

    slice = np.zeros((len(x_range), len(y_range)))

    for i, x in enumerate(x_range):
        for j, y in enumerate(y_range):
            slice[i, j] = image.density_at(x, y, z)

    show_image(slice)

def show_image(slice_array):
    image_array = np.zeros(slice_array.shape, dtype=np.uint8)

    for i in range(image_array.shape[0]):
        for j in range(image_array.shape[1]):
            if slice_array[i, j] > 0.0001:
                image_array[i, j] = 255.255

    img = Image.fromarray(image_array, "L")
    # Display the Numpy array as Image
    img.show()

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
                        help="x coordinate")

    parser.add_argument("-y",
                        "--ysize",
                        type=int,
                        required=True,
                        help="y coordinate")

    parser.add_argument("-z",
                        "--zcoord",
                        type=float,
                        required=True,
                        help="z coordinate")

    return parser.parse_args()

def main():
    """
    run the demo
    """
    args = get_args()
    with mrcfile.open(args.input, mode='r+') as mrc:
        display_cross_section(mli.MRCImage(mrc), args)

if __name__ == "__main__":
    main()
