import argparse
import pathlib
import sys

from copy_mesh import read_tetgen_file
from copy_mesh import remove_suffix

def get_args():
    description = "Verifies that the two meshes provided are identical, ignoring comments."
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("file1",
                        type=pathlib.Path,
                        help="The first input file.")
    
    parser.add_argument("file2",
                        type=pathlib.Path,
                        help="The second input file.")
    
    # Maybe other arguments for mesh optimization later on.
    return parser.parse_args()




def command_validation(args):
    # Verify that the input exists
    if not args.file1.exists():
        return f"Error: file {args.file1} does not exist!"
    
    if not args.file2.exists():
        return f"Error: file {args.file2} does not exist!"
    
    return None

def main():
    args = get_args()
    # Verify arguments
    error_message = command_validation(args)
    if error_message is not None:
        print(error_message, file=sys.stderr)
        return

    # Read files and extract mesh data
    nodes1, faces1, tets1, _  = read_tetgen_file(pathlib.Path(remove_suffix(str(args.file1))))
    nodes2, faces2, tets2, _  = read_tetgen_file(pathlib.Path(remove_suffix(str(args.file2))))

    matching_flag = True

    if nodes1 == nodes2:
        print("Nodes/Points file is matching.")
    else:
        matching_flag = False
        print("Nodes/Point file is not matching.")

    if faces1 == faces2:
        print("Faces file is matching.")
    else:
        matching_flag = False
        print("Faces file is not matching.")

    if tets1 == tets2:
        print("Tets file is matching.")
    else:
        matching_flag = False
        print("Tets file is not matching.")

    if matching_flag:
        print("Final Verdict: Mesh is identical.")
    else:
        print("Final Verdict: Mesh is different.")
    

if __name__ == "__main__":
    main()


    
