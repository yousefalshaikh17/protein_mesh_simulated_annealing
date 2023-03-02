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
import numpy as np

from ffeamesh.vector3 import Vector3

def shortest_sides(node_points, tets, tet_props):
    """
    lists shortest side of each tet
    Args:
        node_points ([NodePoint]): the vertices
        tets ([Tetrahedron4]): the tets
        tet_props (dict(int, [float])): tetgen index => properties of the tet
    """
    for tet in tets.values():
        nodes = []
        nodes.append(node_points[tet.vert0])
        nodes.append(node_points[tet.vert1])
        nodes.append(node_points[tet.vert2])
        nodes.append(node_points[tet.vert3])

        shortest = nodes[0].to_edge_vector(nodes[1]).length()

        for node in nodes[2:]:
            tmp = nodes[0].to_edge_vector(node).length()
            if tmp < shortest:
                shortest = tmp

            tmp = nodes[1].to_edge_vector(node).length()
            if tmp < shortest:
                shortest = tmp

        tmp = nodes[2].to_edge_vector(nodes[3]).length()
        if tmp < shortest:
            shortest = tmp

        tet_props[tet.index].append(shortest)

def tet_area(coords):
    """
    Find tet area by summing the areas of the faces
    Args:
        coords [NodePoint]: the four vertices of the tet
    Returns:
        (float): total surface area of tet
    """
    sides = []
    for coord in coords[1:]:
        sides.append(coords[0].to_edge_array(coord))

    sides.append(coords[1].to_edge_array(coords[2]))
    sides.append(coords[1].to_edge_array(coords[3]))

    area = [np.cross(sides[0], sides[1])]
    area.append(np.cross(sides[1], sides[2]))
    area.append(np.cross(sides[2], sides[0]))
    area.append(np.cross(sides[3], sides[4]))

    area = [np.linalg.norm(x)/2.0 for x in area]

    return sum(area)

def tet_volume(coords):
    """
    Find tet volume by ((side1 x side2).side2)/6
    Args:
        coords [[x, y, z]]: the four vertices of the tet
    Returns:
        (float): the volume of the tet
    """
    sides = []
    for coord in coords[1:]:
        sides.append(coords[0].to_edge_array(coord))

    return abs(np.dot(np.cross(sides[0], sides[1]), sides[2]))/6.0

def areas(node_points, tets, tet_props):
    """
    lists surface area of each tet
    Args:
        node_points ([NodePoint]): the vertices
        tets ([Tetrahedron4]): the tets
        tet_props (dict(int, [float])): tetgen index => properties of the tet
    """
    for tet in tets.values():
        nodes = []
        nodes.append(node_points[tet.vert0])
        nodes.append(node_points[tet.vert1])
        nodes.append(node_points[tet.vert2])
        nodes.append(node_points[tet.vert3])

        tet_props[tet.index].append(tet_area(nodes))

def volumes(node_points, tets, tet_props):
    """
    lists volume of each tet
    Args:
        node_points ([NodePoint]): the vertices
        tets ([Tetrahedron4]): the tets
        tet_props (dict(int, [float])): tetgen index => properties of the tet
    """
    for tet in tets.values():
        nodes = []
        nodes.append(node_points[tet.vert0])
        nodes.append(node_points[tet.vert1])
        nodes.append(node_points[tet.vert2])
        nodes.append(node_points[tet.vert3])

        tet_props[tet.index].append(tet_volume(nodes))

def get_tet_props(node_points, tets):
    """
    lists shortest side of each tet
    Args:
        node_points ([NodePoint]): the vertices
        tets ([Tetrahedron4]): the tets
    Returns:
        tet_props (dict(int, [float])): tetgen index => properties of the tet
    """
    tet_props = {}

    # initialize dict
    for tet in tets.values():
        tet_props[tet.index] = []

    shortest_sides(node_points, tets, tet_props)
    volumes(node_points, tets, tet_props)
    areas(node_points, tets, tet_props)

    return tet_props

def edges_to_area_ratio(nodes):
    """
    find the ratio of the square of the total edge length to the area
    Args:
        nodes (NodePoint): the three vertices
    Returns
        (numpy.float64): L^2/A
    """
    sides = []
    sides.append(nodes[0].to_edge_array(nodes[1]))
    sides.append(nodes[0].to_edge_array(nodes[2]))
    sides.append(nodes[1].to_edge_array(nodes[2]))

    # half the length of the cross product
    area = np.linalg.norm(np.cross(sides[0], sides[1]))/2.0

    length = sum(np.linalg.norm(x) for x in sides)

    return length*length/area

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
        offsets.append(Vector3([tmp_x, tmp_y, tmp_z]))

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

def uniformity(values):
    """
    find the variance of the values devided by the square of the average
    Args:
        values (numbers)
    Return:
        (numpy.float64): variance/average^2
    """
    return np.var(values)/(np.mean(values)**2)