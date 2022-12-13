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
import csv
import pathlib
import collections
import argparse
import sys

## data structure for metadata line of tetgen .node file
## number of points, dimension, number of attributes and boundary markers
NodeMetaData = collections.namedtuple("NodeMetaData", "points, dimension, attributes, bms")

## 3D point with index
NodePoint = collections.namedtuple("NodePoint", "index, x, y, z")

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
                        required=True,
                        help="root name of .node .face & .ele files")

    return parser.parse_args()

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
        meta_data = NodeMetaData(int(row[0]), int(row[1]), int(row[2]), int(row[3]))

        if meta_data.dimension != 3:
            raise ValueError(f"File {input_file} is not 3D.")

        points = []
        for row in reader:
            points.append(NodePoint(int(row[0]), float(row[1]),  float(row[2]), float(row[3]),))

        if len(points) != meta_data.points:
            req = meta_data.points
            act = len(points)
            er_m = f"File {input_file} should have {req} points, but {act} were found!"
            raise ValueError(er_m)

    return meta_data, points

def print_data(meta_data, data):
    """
    print the node data
    """
    print(meta_data)

    for item in data:
        print(item)

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
        meta_data = FaceMetaData(int(row[0]), int(row[1]))

        faces = []
        for row in reader:
            faces.append(Face(int(row[0]), int(row[1]),  int(row[2]), int(row[3]), int(row[4])))

        if len(faces) != meta_data.faces:
            req = meta_data.faces
            act = len(faces)
            er_m = f"File {input_file} should have {req} points, but {act} were found!"
            raise ValueError(er_m)

    return meta_data, faces

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
        meta_data = TetMetaData(int(row[0]), int(row[1]), int(row[2]))

        if meta_data.nodes != 4:
            raise ValueError("This reader can only handle four nodes per tetrahedron.")

        tets = []
        for row in reader:
            if meta_data.ra == 0:
                tets.append(Tetrahedron4(int(row[0]),
                                         int(row[1]),
                                         int(row[2]),
                                         int(row[3]),
                                         int(row[4]),
                                         None))
            else:
                tets.append(Tetrahedron4(int(row[0]),
                                             int(row[1]),
                                             int(row[2]),
                                             int(row[3]),
                                             int(row[4],
                                             int(row[5]))))

        if len(tets) != meta_data.tets:
            req = meta_data.tets
            act = len(tets)
            er_m = f"File {input_file} should have {req} points, but {act} were found!"
            raise ValueError(er_m)

    return meta_data, tets

def main():
    """
    run the program
    """
    args = get_args()
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

    try:
        node_data, node_points = read_node_file(node_file)
        print_data(node_data, node_points)

        face_data, faces = read_face_file(face_file)
        print_data(face_data, faces)

        tets_data, tets = read_tet_file(tets_file)
        print_data(tets_data, tets)
    except ValueError as error:
        print(error, file=sys.stderr)
        return

if __name__ == "__main__":
    main()
