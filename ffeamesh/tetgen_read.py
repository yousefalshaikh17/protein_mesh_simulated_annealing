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

 @author: jonathan pickering, joanna leng 20 Dec 22
"""
import csv
import collections
import operator
import sys
import numpy as np

## data structure for metadata line of tetgen .node file
## number of points, dimension, number of attributes and boundary markers
NodeMetaData = collections.namedtuple("NodeMetaData", "points, dimension, attributes, bms")

## 3D point with index
_NodePoint = collections.namedtuple("_NodePoint", "index, x, y, z")

class NodePoint(_NodePoint):
    """
    a point with an index
    """

    def to_edge(self, rhs):
        """
        vector from self to rhs
        Args:
            rhs (NodePoint): the 'to' vector
        Returns:
            EdgeVector
        """
        return StlVector([rhs.x - self.x, rhs.y - self.y, rhs.z - self.z])

    def to_stl_vertex(self):
        """
        to text stl vertex format
        """
        return f"vertex {self.x:e} {self.y:e} {self.z:e}"

class StlVector(list):
    """
    a vector represnting an edge or normal
    """

    def dot(self, rhs):
        """
        dot (inner) product self.rhs
        Args:
            rhs (NodePoint): right hand side of dot
        Returns
            float
        """
        return self[0]*rhs[0] + self[1]*rhs[1] + self[2]*rhs[2]

    def cross(self, rhs):
        """
        cross product self x rhs
        Args:
            rhs (NodePoint): right hand side of cross
        Returns:
            (float, float, float)
        """
        return StlVector(np.cross(self, rhs))

    def normalize(self):
        """
        make a normalized copy of objects vector
        Returns
            [float]
        """
        length = np.linalg.norm(self)
        tmp = np.divide(self, length)
        return StlVector(tmp)

    def isclose(self, rhs):
        """
        floating point comparison of two vectors
        Args:
            rhs (StlVector): comaprison object
        Returns:
            True if all comonants pass numpy isclose
        """
        results = np.isclose(self, rhs)
        return all(x for x in results)

    def to_stl_normal(self):
        """
        make a stl fromat normal vector
        """
        return f"normal {self[0]:e} {self[1]:e} {self[2]:e}"


def make_stl_normal(node0, node1, node2):
    """
    make an stl normal from nodes assume vertices are
    listed in counter-clock-wise order from outside
    """
    vec0 = node0.to_edge(node1)
    vec1 = node1.to_edge(node2)

    return vec0.cross(vec1).normalize()

## data structure for metadata line of tetgen .node file
## number of faces and boundary markers
FaceMetaData = collections.namedtuple("FaceMetaData", "faces, bm")

## a face
Face = collections.namedtuple("Face", "index, vert0, vert1, vert2, bm")

## data structure for metadata line of tetgen .node file
## number of tets, number of nodes pre tet and region attribute
TetMetaData = collections.namedtuple("TetMetaData", "tets, nodes, ra")

## a face
Tetrahedron4 = collections.namedtuple("Tetrahedron4", "index, vert0, vert1, vert2, vert3, ra")

def make_decomment(comment):
    """
    make a decomment function
    Args:
        comment (str): the comment string
    """

    if comment is None or comment.isspace() or comment == '':
        return None

    def decomment(csvfile):
        """
        make iterator to remove lines starting with the comment symbols
        Args:
            csvfile ()
        Returns:
            generator returning non-comment file rows as (str)
        """
        for row in csvfile:
            if len(row)>0 and not row.strip().startswith(comment):
                row = row[:row.find(comment)].strip()
                if len(row)>0:
                    yield row

    return decomment

def read_node_file(input_file):
    """
    read a tetgen nodes file
    Args:
        input_file (pathlib.Path): source
    Return
        (NodeMetaData): the file's metadata
        ([Point]): the points
    Throws:
        ValueError if problem
    """
    decomment = make_decomment("#")

    with input_file.open('r', newline='') as file:
        reader = csv.reader(decomment(file), delimiter=' ')
        row = next(reader)
        row = [x for x in row if x != '']
        meta_data = NodeMetaData(int(row[0]), int(row[1]), int(row[2]), int(row[3]))

        if meta_data.dimension != 3:
            raise ValueError(f"File {input_file} is not 3D.")

        points = {}
        for row in reader:
            row = [x for x in row if x != '']
            index = int(row[0])
            points[index] = NodePoint(index, float(row[1]),  float(row[2]), float(row[3]))

        if len(points) != meta_data.points:
            req = meta_data.points
            act = len(points)
            er_m = f"File {input_file} should have {req} points, but {act} were found!"
            raise ValueError(er_m)

    return meta_data, points

def read_face_file(input_file):
    """
    read a tetgen faces file
    Args:
        input_file (pathlib.Path): source
    Return
        (FaceMetaData): the file's metadata
        ([Face]): the faces
    Throws:
        ValueError if problem
    """
    decomment = make_decomment("#")

    with input_file.open('r', newline='') as file:
        reader = csv.reader(decomment(file), delimiter=' ')
        row = next(reader)
        row = [x for x in row if x != '']
        meta_data = FaceMetaData(int(row[0]), int(row[1]))

        faces = []
        for row in reader:
            row = [x for x in row if x != '']
            faces.append(Face(int(row[0]), int(row[1]),  int(row[2]), int(row[3]), int(row[4])))

        if len(faces) != meta_data.faces:
            req = meta_data.faces
            act = len(faces)
            er_m = f"File {input_file} should have {req} points, but {act} were found!"
            raise ValueError(er_m)

    return meta_data, sorted(faces, key=operator.attrgetter('index'))

def read_tet_file(input_file):
    """
    read a tetgen ele file
    Args:
        input_file (pathlib.Path): source
    Return
        (TetMetaData): the file's metadata
        ([Tetrahedron4]): the tets
    Throws:
        ValueError if problem
    """
    decomment = make_decomment("#")

    with input_file.open('r', newline='') as file:
        reader = csv.reader(decomment(file), delimiter=' ')
        row = next(reader)
        row = [x for x in row if x != '']
        meta_data = TetMetaData(int(row[0]), int(row[1]), int(row[2]))

        if meta_data.nodes != 4:
            raise ValueError("This reader can only handle four nodes per tetrahedron.")

        tets = {}
        for row in reader:
            row = [x for x in row if x != '']
            index = int(row[0])
            if meta_data.ra == 0:
                tets[index] = Tetrahedron4(index,
                                           int(row[1]),
                                           int(row[2]),
                                           int(row[3]),
                                           int(row[4]),
                                           None)
            else:
                tets[index] = Tetrahedron4(index,
                                           int(row[1]),
                                           int(row[2]),
                                           int(row[3]),
                                           int(row[4],
                                           int(row[5])))

        if len(tets) != meta_data.tets:
            req = meta_data.tets
            act = len(tets)
            er_m = f"File {input_file} should have {req} points, but {act} were found!"
            raise ValueError(er_m)

    return meta_data, tets

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
        print(f"endsolid", file=out_file)

def print_stl_facet(node0, node1, node2, out_stream=sys.stdout):
    """
    convert three nodes to an stl facet and print
    Args:
        node1 (NodePoint)
        node1 (NodePoint)
        node1 (NodePoint)
        out_stream (iotextstream): output destination
    """
    normal = make_stl_normal(node0, node1, node2)

    print(f" facet {normal.to_stl_normal()}", file=out_stream)
    print(" outer loop", file=out_stream)
    print(f"   {node0.to_stl_vertex()}", file=out_stream)
    print(f"   {node1.to_stl_vertex()}", file=out_stream)
    print(f"   {node2.to_stl_vertex()}", file=out_stream)
    print(" endloop", file=out_stream)
    print(" endfacet", file=out_stream)

def measure_excentricity(node_points, tets):
    """
    build tets and check excentricity
    Args:
        node_points ([NodePoint]): the vertices
        tets ([Tetrahedron4]): the tets
    Returns:
        ([int])
    """
    ratios = {}
    for tet in tets.values():
        nodes = []
        nodes.append(node_points[tet.vert0])
        nodes.append(node_points[tet.vert1])
        nodes.append(node_points[tet.vert2])
        nodes.append(node_points[tet.vert3])

        inertia_tensor = make_inertia_tensor(nodes)
        e_vals, _ = np.linalg.eig(inertia_tensor)

        scales = [abs(x) for x in e_vals]
        scales.sort()

        if scales[0] > 0.0:
            ratios[tet.index] = scales[2]/scales[0]
        else:
            ratios[tet.index] = sys.float_info.max

    return ratios

def make_inertia_tensor(nodes):
    """
    construct the inertia tensor of a set of nodes,
    (assume unit mass each node point unit mass)
    Args:
        nodes ([NodePoint])
    Returns:
        (numpy.array 3x3 float)
    """
    ctr_x, ctr_y, ctr_z = find_centre(nodes)

    offsets = []
    for node in nodes:
        tmp_x = node.x - ctr_x
        tmp_y = node.y - ctr_y
        tmp_z = node.z - ctr_z
        offsets.append(StlVector([tmp_x, tmp_y, tmp_z]))

    inertia_tensor = np.zeros([3, 3], dtype=float)
    for vec in offsets:
        # leading diagonal
        squares = [x*x for x in vec]
        inertia_tensor[0, 0] += squares[1] + squares[2]
        inertia_tensor[1, 1] += squares[0] + squares[2]
        inertia_tensor[1, 1] += squares[0] + squares[1]

        # uppar triangular
        inertia_tensor[0, 1] -= vec[0]*vec[1]
        inertia_tensor[0, 2] -= vec[0]*vec[2]
        inertia_tensor[1, 2] -= vec[1]*vec[2]

    # lower triangular (symmetric with uppar)
    inertia_tensor[1, 0] = inertia_tensor[0, 1]
    inertia_tensor[2, 0] = inertia_tensor[0, 2]
    inertia_tensor[2, 1] = inertia_tensor[1, 2]

    return inertia_tensor

def find_centre(nodes):
    """
    find the centre of a set of nodes
    """
    ctr_x = 0.0
    ctr_y = 0.0
    ctr_z = 0.0

    for node in nodes:
        ctr_x += node.x
        ctr_y += node.y
        ctr_z += node.z

    length = float(len(nodes))

    return ctr_x/length, ctr_y/length, ctr_z/length

#########################################################################