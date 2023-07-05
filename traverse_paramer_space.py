"""
Inputs:
    Prob mutate 0-1
    max mutate size 0-2.5
    all/ single axis?
    weights 0-inf
    cooling type, speed

Outputs:
    Shape (average, min, max)
    Size (average, min, max)
    isosurface match {make image using paraview}
"""
import argparse
import os
import pathlib

def get_args():
    """
    get the command line arguments
    Returns:
        (argparse.namespace)
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("-f",
                        "--params_files",
                        type=pathlib.Path,
                        required=True,
                        help="file holding list of sim_anneal parameters files")

    return parser.parse_args()


def main():
    args = get_args()
    print(args)

    s = "python .\\ffeamesh\\scripts\\sim_anneal.py file "
    s += f" -f {args.params_files} "
    print(s)
    os.system(s)
    print("\tfin")

if __name__ == "__main__":
    main()
