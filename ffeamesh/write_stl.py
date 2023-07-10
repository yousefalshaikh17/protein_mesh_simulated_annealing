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

 @author: jonathan pickering, joanna leng 08 Feb 23
"""
# set up linting
# pylint: disable = import-error
import pathlib
import sys

def to_stl_normal_string(vec):
    """
    make a stl fromat normal vector
    Args:
        vec (Vector3)
    """
    return f"normal {vec[0]:e} {vec[1]:e} {vec[2]:e}"

def make_vector_normal(node0, node1, node2):
    """
    make an stl normal from nodes assume vertices are
    listed in counter-clock-wise order from outside
    """
    vec0 = node0.to_edge(node1)
    vec1 = node1.to_edge(node2)

    return vec0.cross(vec1).normalize()

def node_to_stl_vertex_string(node):
    """
    to text stl vertex format
    """
    return f"vertex {node.x:e} {node.y:e} {node.z:e}"

def write_stl(nodes, faces, stl_name):
    """
    write surface to STL file
    Args:
        name_root (pathlib.Path): root name of system
        nodes ([NodePoint])
        faces ([Face])
    """
    if stl_name.suffix != '.stl':
        stl_name = pathlib.Path(str(stl_name)+".stl")

    print_stl_file(nodes, faces, stl_name)
    print(f"to_stl: {stl_name}, {len(faces)} faces on {len(nodes)} nodes")

def print_stl_file(nodes, faces, stl_name):
    """
    write an stl file
    """
    with stl_name.open('w') as out_file:
        print(f"solid {stl_name.stem}", file=out_file)
        for face in faces:
            print_stl_facet(nodes[face.vert0], nodes[face.vert1], nodes[face.vert2], out_file)
        print("endsolid", file=out_file)

def print_stl_facet(node0, node1, node2, out_stream=sys.stdout):
    """
    convert three nodes to an stl facet and print
    Args:
        node1 (NodePoint)
        node1 (NodePoint)
        node1 (NodePoint)
        out_stream (iotextstream): output destination
    """
    normal = make_vector_normal(node0, node1, node2)

    print(f" facet {to_stl_normal_string(normal)}", file=out_stream)
    print(" outer loop", file=out_stream)
    print(f"   {node0.to_stl_vertex_string()}", file=out_stream)
    print(f"   {node1.to_stl_vertex_string()}", file=out_stream)
    print(f"   {node2.to_stl_vertex_string()}", file=out_stream)
    print(" endloop", file=out_stream)
    print(" endfacet", file=out_stream)
