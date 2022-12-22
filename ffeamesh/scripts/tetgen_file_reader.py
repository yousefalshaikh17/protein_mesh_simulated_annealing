"""
 This file is part of the FFEA simulation package

 Copyright (c) by the Theory and Development FFEA teams,
 as they appear in the README.md file.

 FFEA is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 FFEA is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with FFEA.  If not, see <http://www.gnu.org/licenses/>.

 To help us fund FFEA development, we humbly ask that you cite
 the research papers on the package.

 @author: jonathan pickering, joanna leng 13 Dec 22
"""
import pathlib
import argparse
import sys

import ffeamesh.tetgen_read as tr

def make_test_nodes():
    """
    make some test nodes
    """
    node0 = tr.NodePoint(0, 10.867435455322266, 272.447998046875, -69.545654296875)
    node1 = tr.NodePoint(0, 10.867435455322266, 257.23294067382813, -56.546279907226563)
    node2 = tr.NodePoint(0, -6.920966625213623, 272.447998046875, -56.546279907226563)

    return node0, node1, node2

def demo_print_stl():
    """
    test the construction of edge vectors & their cross product
    Args:

    """
    node0, node1, node2 = make_test_nodes()
    tr.print_stl_facet(node0, node1, node2)



def test_vector():
    """
    test the construction of edge vectors & their cross product
    """
    node0, node1, node2 = make_test_nodes()
    normal = tr.make_stl_normal(node0, node1, node2)

    if normal.isclose(tr.StlVector([-0.48567734492792791, -0.56782065825723671, -0.6646030519641607])):
        print("Test vector: Pass")
    else:
        print("Test vector: Fail")

#TODO move to tests
def make_test_vectors():
    """
    make a set of test nodes
    returns
        NodePoint, NodePoint, NodePoint, NodePoint
    """
    t_x = tr.StlVector([1.0, 0.0, 0.0])
    t_y = tr.StlVector([0.0, 1.0, 0.0])
    t_z = tr.StlVector([0.0, 0.0, 1.0])
    t_t = tr.StlVector([2.0, 3.0, 4.0])

    return t_x, t_y, t_z, t_t

def test_cross():
    """
    test the cross product
    """
    t_x, _, t_z, t_t = make_test_vectors()

    r_x, r_y, r_z = t_x.cross(t_z)

    if r_x != 0.0 and r_y != 1.0 and r_z != 0.0:
        print(f"Cross product: Fail {t_x}x{t_z} = ({r_x}, {r_y}, {r_z})")
        return
    else:
        print("Cross product: Pass")

    r_x, r_y, r_z = t_z.cross(t_x)

    if r_x != 0.0 and r_y != -1.0 and r_z != 0.0:
        print(f"Cross product: Fail {t_z}x{t_x} = ({r_x}, {r_y}, {r_z})")
        return
    else:
        print("Cross product: Pass")

    r_x, r_y, r_z = t_x.cross(t_t)

    if r_x != 0.0 and r_y != -4.0 and r_z != 3.0:
        print(f"Cross product: Fail {t_x}x{t_t} = ({r_x}, {r_y}, {r_z})")
        return
    else:
        print("Cross product: Pass")

def test_dot():
    """
    test the dot product
    """
    t_x, t_y, t_z, t_t = make_test_vectors()

    if t_x.dot(t_y) != 0.0:
        print("Dot product: Fail {t_x}.{t_y} != 0.0")
        return
    else:
        print("Dot product: Pass")

    if t_x.dot(t_z) != 0.0:
        print("Dot product: Fail {t_x}.{t_z} != 0.0")
        return
    else:
        print("Dot product: Pass")

    if t_y.dot(t_z) != 0.0:
        print("Dot product: Fail {t_y}.{t_z} != 0.0")
        return
    else:
        print("Dot product: Pass")

    if t_x.dot(t_t) != 2.0:
        print("Dot product: Fail {t_x}.{t_t} != 2.0")
        return
    else:
        print("Dot product: Pass")

    if t_y.dot(t_t) != 3.0:
        print("Dot product: Fail {t_y}.{t_t} != 3.0")
        return
    else:
        print("Dot product: Pass")

    if t_z.dot(t_t) != 4.0:
        print("Dot product: Fail {t_z}.{t_t} != 4.0")
        return
    else:
        print("Dot product: Pass")

def print_data(root_name, nodes_data, faces_data, tets_data):
    """
    print the node data
    """
    out_s = f"System {root_name}: "
    out_s += f"{nodes_data.points} vertices, "
    out_s += f"{tets_data.tets} tetrahedra, "
    out_s += f"{faces_data.faces} faces"
    print(out_s)

def get_args():
    """
    get command line arguments
    Returns
        (argparse.namespace)
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("-r",
                        "--root_name",
                        type=pathlib.Path,
                        required=False,
                        help="root name of .node .face & .ele files")

    parser.add_argument("-s",
                        "--stl_name",
                        type=pathlib.Path,
                        required=False,
                        help="file for stl outut if required")

    parser.add_argument("-e",
                        "--tets_excentricity",
                        action="store_true",
                        help="if set compute the excentricity of the tets")

    parser.add_argument("-t",
                        "--test",
                        action="store_true",
                        help="if set run the test harness")

    return parser.parse_args()

def process_files(node_file, face_file, tets_file, args):
    """
    process the three files
    """
    try:
        nodes_data, node_points = tr.read_node_file(node_file)
        faces_data, faces       = tr.read_face_file(face_file)
        tets_data, tets         = tr.read_tet_file(tets_file)

        print_data(args.root_name, nodes_data, faces_data, tets_data)

        if args.stl_name is not None:
            tr.write_stl(node_points, faces, args.stl_name)


        if args.tets_excentricity:
            ratios = tr.measure_excentricity(node_points, tets)
            maximum = max(ratios.values())
            print(f"Max: {maximum}")#, tet number {max_index+1}")

    except ValueError as error:
        print(error, file=sys.stderr)
        return

def main():
    """
    run the program
    """
    args = get_args()

    if args.test:
        test_dot()
        test_cross()
        test_vector()
        demo_print_stl()
        return

    if args.root_name is None:
        return

    node_file = args.root_name.with_suffix(".1.node")
    face_file = args.root_name.with_suffix(".1.face")
    tets_file = args.root_name.with_suffix(".1.ele")

    if not node_file.exists():
        print(f"File {node_file} doesn't exist", file=sys.stderr)
        return

    if not face_file.exists():
        print(f"File {face_file} doesn't exist", file=sys.stderr)
        return

    if not tets_file.exists():
        print(f"File {tets_file} doesn't exist", file=sys.stderr)
        return

    process_files(node_file, face_file, tets_file, args)

if __name__ == "__main__":
    main()
