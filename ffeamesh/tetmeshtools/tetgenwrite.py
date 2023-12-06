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

Authors: Joanna Leng, Jonathan Pickering, - University of Leeds (UK)
"""

def write_tetgen_elements(output_file, tets, comment=""):
    """
    write ffea elements file .ele
    First line: <# of tetrahedra> <nodes per tetrahedron> <# of attributes>
    Remaining lines list of # of tetrahedra:
    <tetrahedron #> <node> <node> <node> <node> ... [attributes]
    Args:
        output_file (pathlib.Path): root name of file, will have .node added
        tets ({index => Tetrahedron4}):
        comment (str): any user comment
    """
    with open(output_file.with_suffix(".1.ele"), "w", encoding="utf-8") as ele:
        ele.write(f"{len(tets)} 4 0\n")
        for index in tets:
            tet = tets[index]
            ele.write(f"{index} {tet.vert0} {tet.vert1} {tet.vert2} {tet.vert3}\n")
        ele.write(f"# {comment}")

def write_tetgen_nodes(output_file, points, comment=""):
    """
    write ffea .node file
    First line: <num points> <dimension (3)> <num attributes> <num boundary markers (0 or 1)>
    Remaining lines list # of points:
    <point #> <x> <y> <z>
    Args:
        output_file (pathlib.Path): root name of file, will have .node added
        points ({index => NodePoint}):
        comment (str): any user comment
    """
    with open(output_file.with_suffix(".1.node"), "w", encoding="utf-8") as node:
        node.write(f"{len(points)} 3 0 0\n")
        for index in points:
            point = points[index]
            node.write(f"{index} {point.x} {point.y} {point.z}\n")
        node.write(f"# {comment}")

def write_tetgen_faces(output_file, faces, comment=""):
    """
    write .face file
    First line: <# of faces> <boundary marker (0 or 1)>
    Remaining lines list of # of faces:
    <face #> <node> <node> <node> [boundary marker]
    Args:
        output_file (pathlib.Path): root name of file, will have .face added
        faces ({index => Face}): connectivity of face polygons
        comment (str): any user comment
    """
    with open(output_file.with_suffix(".1.face"), "w", encoding="utf-8") as face_file:
        face_file.write(f"{len(faces)} 1\n")
        for index in faces:
            face = faces[index]
            face_file.write(f"{index} {face.vert0} {face.vert1} {face.vert2} -1\n")
        face_file.write(f"# {comment}")
