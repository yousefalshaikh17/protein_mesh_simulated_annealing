import argparse
import pathlib
import sys
import tetmeshtools.tetprops as tp
import tetmeshtools.meshtools.tetgenread as tr
import tetmeshtools.meshtools.tetgenwrite as tw


def get_args():
    description = "Makes a copy of the input tetmesh and saves it as the output name."
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("-i",
                        "--input",
                        type=pathlib.Path,
                        required=True,
                        help="input file")
    
    parser.add_argument("-o",
                        "--output",
                        type=pathlib.Path,
                        required=True,
                        help="output file name root, (no suffix)")
    
    # Maybe other arguments for mesh optimization later on.
    return parser.parse_args()


def remove_suffix(fname):
    for suffix in [".1.node", ".1.ele", ".1.face"]:
            if fname.endswith(suffix):
                return fname.replace(suffix, "", 1)

def read_tetgen_file(root_name):
    tet_files = tr.make_and_test_tetgen_files(root_name)
    try:
            _, nodes = tr.read_node_file(tet_files[0])
            _, faces = tr.read_face_file(tet_files[1])
            _, tets  = tr.read_tet_file(tet_files[2])
            tet_props = tp.get_tet_props(nodes, tets)
    except ValueError as error:
            print(error)
            return
    return nodes, faces, tets, tet_props
    
def write_tetgen_file(root_name, nodes, faces, tets, comment=None):
    tw.write_tetgen_faces(root_name, faces, comment)
    tw.write_tetgen_elements(root_name, tets, comment)
    tw.write_tetgen_nodes(root_name, nodes, comment)

def print_mesh_data(nodes, faces, tet_props):
    print("Nodes:")
    print(f"{'Index':<10}{'X':<15}{'Y':<15}{'Z':<15}")
    print('-' * 50)
    for i,node in nodes.items():
        if i > 5:
             break
        print(f"{node.index:<10}{node.x:<15.5f}{node.y:<15.5f}{node.z:<15.5f}")

    print("\n\n")

    print("Faces:")
    print(f"{'Index':<10}{'Vertix0':<15}{'Vertix1':<15}{'Vertix2':<15}{'Boundary Markers':<15}")
    print('-' * 75)
    for i,face in faces.items():
        if i > 5:
             break
        print(f"{face.index:<10}{face.vert0:<15}{face.vert1:<15}{face.vert2:<15}{face.bm:<15.5f}")

    print("\n\n")

    # print("Elements:")
    # print(f"{'Index':<10}{'Vertix0':<15}{'Vertix1':<15}{'Vertix2':<15}{'Vertix3':<15}{'RA':<15}")
    # print('-' * 75)
    # for i in tets:
    #     tet = tets[i]
    #     print(f"{tet.index:<10}{tet.vert0:<15.5f}{tet.vert1:<15.5f}{tet.vert2:<15.5f}{tet.vert3:<15.5f}" f"{f'(tet.ra:<15.5f)' if tet.ra is not None else 'N/A'}")

    print("Elements:")
    print(f"{'TetGenIndex':<13}{'Shortest Side':<15}{'Volume':<15}{'Surface Area':<15}{'Shape Factor':<15}")
    print('-' * 75)
    for i, tet_prop in tet_props.items():
        if i > 5:
             break
        print(f"{tet_prop[0]:<13}{tet_prop[1]:<15}{tet_prop[2]:<15.5f}{tet_prop[3]:<15.5f}{tet_prop[4]:<15.5f}")


def command_validation(args):
    """
    Validates commands.
    """
    # Verify that the input exists
    if not args.input.exists():
        return f"Error: file {args.input} does not exist!"
    return None

# Unit testing for surface tetrahedron
def verify_surface_tetrahedron(surface_tets, tets, faces):
    # Number of surface tetrahedron should match number of faces.
    if len(faces) != len(surface_tets):
         return False
    # Verify that the face coordinates are also included in tetrahedron coordinates. 
    for face in faces.values():
        # Verify that the face has a key in the surface tetrahedron dictionary
        if surface_tets[face.index] is None:
             return False
        # Verify that the face coordinates are also included in tetrahedron coordinates. 
        tet = tets[surface_tets[face.index]]
        if not {face.vert0, face.vert1, face.vert2}.issubset({tet.vert0, tet.vert1, tet.vert2, tet.vert3}):
            return False
    return True

import time # remove later. *****
def main():
    args = get_args()
    # Verify arguments
    error_message = command_validation(args)
    if error_message is not None:
        print(error_message, file=sys.stderr)
        return

    # Read file and extract mesh data
    nodes, faces, tets, _ = read_tetgen_file(pathlib.Path(remove_suffix(str(args.input))))

    # Print mesh
    print_mesh_data(nodes, faces, tets)

    # TODO: Process the mesh
    # strt_time = time.time()
    # surface_tets = get_surface_tets(tets, faces)
    # print(time.time()-strt_time)
    

    # print(len(surface_tets))
    # print(len(faces))

    # Verify surface tets (*** remove)
    # if verify_surface_tetrahedron(surface_tets, tets, faces):
    #     print("Surface tetrahedrons verified.")
    # else:
    #      print("Surface tetrahedrons verification failed.")

    # Save mesh data as a file
    write_tetgen_file(pathlib.Path(str(args.output)), nodes, faces, tets)
    print("\n\nSaved mesh files.")


if __name__ == "__main__":
    main()