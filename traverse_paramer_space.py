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

def run_file(command):
    """
    run a single file
    Args:
        source_file (str): the path to the parameters file
    """
    print(f"run:\t{command}")
    os.system(command)

def main():
    """
    run the script
    """
    args = get_args()
    with args.params_files.open('r') as file:
        lines = file.readlines()
        for line in lines:
            run_file(line.strip())

    print("\t>>END<<")

if __name__ == "__main__":
    main()
