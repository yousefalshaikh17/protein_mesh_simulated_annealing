import argparse
import sys
import mrcfile
import pathlib
import tetmeshtools.mrclininterpolate as mi
from tetmeshtools.scripts.simulated_annealing.copy_mesh import read_tetgen_file, remove_suffix

def get_args():
    description = "Calculates the density of each node and face given an MRC file and a tetmesh file."
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("-m",
                        "--mrcfile",
                        type=pathlib.Path,
                        required=True,
                        help="input MRC file")
    
    parser.add_argument("-t",
                        "--tetmeshfile",
                        type=pathlib.Path,
                        required=True,
                        help="input tetmesh file")
    # Maybe other arguments for mesh optimization later on.
    return parser.parse_args()


def command_validation(args):
    # Verify that the inputs exist
    if not args.mrcfile.exists():
        return f"Error: file {args.mrcfile} does not exist!"
    if not args.tetmeshfile.exists():
        return f"Error: file {args.tetmeshfile} does not exist!"
    return None

def get_density_at_point(image, x,y,z):
    coords = image.to_coords(x,y,z)
    if not isinstance(coords, int):
        try:
            return image.linear_interp(coords)
        except Exception as E:
            print(f"Error calculating density: {E}")
    return None

def compute_node_densities(image, nodes):
    node_densities = dict()
    for node in nodes.values():
        node_densities[node.index] = get_density_at_point(image, node.x, node.y, node.z)
    return node_densities

def compute_face_densities(faces, node_densities):
    face_densities = dict()
    for face in faces.values():
        # Calculate average between all three nodes in the faces
        face_densities[face.index] = (node_densities[face.vert0] + node_densities[face.vert1] + node_densities[face.vert2]) / 3
    return face_densities

def compute_element_densities(tets, node_densities):
    element_densities = dict()
    for tet in tets.values():
        # Calculate average between all four nodes in the tetrahedron
        element_densities[tet.index] = (node_densities[tet.vert0] + node_densities[tet.vert1] + node_densities[tet.vert2] + node_densities[tet.vert3]) / 4
    return element_densities

def print_node_densities(nodes, node_densities):
    print("\n\n")
    print("Node Densities:\n")
    print(f"{'Index':<10}{'X':<15}{'Y':<15}{'Z':<15}{'Density':<15}")
    print('-' * 70)
    for i, node in nodes.items():
        if i > 5:
            break
        print(f"{node.index:<10}{node.x:<15.5f}{node.y:<15.5f}{node.z:<15.5f}{node_densities[i]:<15.8f}")

def print_face_densities(faces, face_densities):
    print("\n\n")
    print("Face Densities:\n")
    print(f"{'Index':<10}{'Vertix0':<10}{'Vertix1':<10}{'Vertix2':<10}{'Density':<15}")
    print('-' * 55)
    for i, face in faces.items():
        if i > 5:
            break
        print(f"{face.index:<10}{face.vert0:<10}{face.vert1:<10}{face.vert2:<10}{face_densities[i]:<15.8f}")

def print_element_densities(tets, tet_densities):
    print("\n\n")
    print("Tetrahedron Densities:\n")
    print(f"{'Index':<10}{'Vertix0':<10}{'Vertix1':<10}{'Vertix2':<10}{'Vertix3':<10}{'Density':<15}")
    print('-' * 65)
    for i, tet in tets.items():
        if i > 5:
            break
        print(f"{tet.index:<10}{tet.vert0:<10}{tet.vert1:<10}{tet.vert2:<10}{tet.vert3:<10}{tet_densities[i]:<15.8f}")

def main():
    # Verify arguments
    args = get_args()
    error_message = command_validation(args)
    if error_message is not None:
        print(error_message, file=sys.stderr)
        return
    
    # load tetgen mesh 
    nodes, faces, tets, _ = read_tetgen_file(pathlib.Path(remove_suffix(str(args.tetmeshfile))))
    path = pathlib.Path(str(args.mrcfile))

    with mrcfile.mmap(path, mode='r+') as mrc:
        mrc_image = mi.MRCImage(mrc)
    
    # Iterate through all nodes and compute/print densities for them
    node_densities = compute_node_densities(mrc_image, nodes)
    print_node_densities(nodes, node_densities)

    # Use previous node density dictionary to calculate face densities and print them.
    # Face density is calculated as the average density across all three nodes.
    # Use of a dictionary makes compute more efficient
    face_densities = compute_face_densities(faces, node_densities)
    print_face_densities(faces, face_densities)

    # Calculate tetrahedron densities and print them to console.
    element_densities = compute_element_densities(tets, node_densities)
    print_element_densities(tets, element_densities)

if __name__ == "__main__":
    main()