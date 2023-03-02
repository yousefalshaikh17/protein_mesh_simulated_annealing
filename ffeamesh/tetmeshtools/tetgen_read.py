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

@copyright 2023
@author: j.h.pickering@leeds.ac.uk and j.leng@leeds.ac.uk
"""
import csv
import collections
import operator

from ffeamesh.vector3 import Vector3

## data structure for metadata line of tetgen .node file
## number of points, dimension, number of attributes and boundary markers
NodeMetaData = collections.namedtuple("NodeMetaData", "points, dimension, attributes, bms")

## 3D point with index base
_NodePoint = collections.namedtuple("_NodePoint", "index, x, y, z")

class NodePoint(_NodePoint):
    """
    extend 3D point
    """

    def to_edge_array(self, rhs):
        """
        vector from self to rhs
        Args:
            rhs (NodePoint): the 'to' vector
        Returns:
            EdgeVector
        """
        return [rhs.x - self.x, rhs.y - self.y, rhs.z - self.z]

    def to_edge_vector(self, rhs):
        """
        vector from self to rhs
        Args:
            rhs (NodePoint): the 'to' vector
        Returns:
            EdgeVector
        """
        return Vector3(self.to_edge_array(rhs))

## data structure for metadata line of tetgen .node file
## number of faces and boundary markers
FaceMetaData = collections.namedtuple("FaceMetaData", "faces, bm")

## a face
Face = collections.namedtuple("Face", "index, vert0, vert1, vert2, bm")

## data structure for metadata line of tetgen .node file
## number of tets, number of nodes pre tet and region attribute
TetMetaData = collections.namedtuple("TetMetaData", "tets, nodes, ra")

## a tetrahedron
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
        ({NodePoint}): the points
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
