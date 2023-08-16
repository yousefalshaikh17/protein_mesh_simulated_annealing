import argparse
import pathlib

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
    print(args)

if __name__ == "__main__":
    main()
