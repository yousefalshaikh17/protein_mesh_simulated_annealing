import argparse
import pathlib
import mrcfile

def get_args():
    """
    get the command line arguments
        Returns
            (argparse.namespace)
    """
    parser = argparse.ArgumentParser(description="""make simple mrc image files
        for use in testing""")

    parser.add_argument("-i",
                        "--input",
                        type=pathlib.Path,
                        required=True,
                        help="input file")

    return parser.parse_args()

def main():
    """
    run the script
    """
    args = get_args()

    if not args.input.exists():
        print(f"Error file {args.input} does not exist.")
        return

    with mrcfile.mmap(args.input, 'r') as mrc:
        for item in mrc.header.dtype.names:
            print(f"{item} => {mrc.header[item]}")

if __name__ == "__main__":
    main()
