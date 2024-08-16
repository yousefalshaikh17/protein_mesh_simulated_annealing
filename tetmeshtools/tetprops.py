"""
functions for finding the geometric properties of tetrahedrons.

---------------------------------------

Licensed under the GNU General Public License, Version 3.0 (the "License"); you
may not use this file except in compliance with the License. You may obtain a
copy of the License at <https://www.gnu.org/licenses/gpl-3.0.html>.

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

This work was funded by Joanna Leng's EPSRC funded RSE Fellowship (EP/R025819/1)

@copyright 2020
@author: j.h.pickering@leeds.ac.uk and j.leng@leeds.ac.uk
"""
# set up linting conditions
# pylint: disable = import-error
import numpy as np
import tetmeshtools.vector3 as v3

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

    # length = area of trapesium = area triangle * 2
    tmp = np.cross(sides[0], sides[1])
    return abs(np.dot(tmp, sides[2]))/6.0

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

def total_volume(node_points, tets):
    """
    total volume of tests
    Args:
        node_points ([NodePoint]): the vertices
        tets ([Tetrahedron4]): the tets
    """
    total = 0.0
    for tet in tets.values():
        nodes = []
        nodes.append(node_points[tet.vert0])
        nodes.append(node_points[tet.vert1])
        nodes.append(node_points[tet.vert2])
        nodes.append(node_points[tet.vert3])

        total += tet_volume(nodes)

    return total

def get_tet_props(node_points, tets):
    """
    gets sortest side, volumes & areas
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

def area_of_triangle(nodes, signed=True):
    """
    find the area of a triangle
    Args:
        nodes (NodePoint): the three vertices
    Returns
        (numpy.float64)
    """
    sides = []
    sides.append(nodes[0].to_edge_array(nodes[1]))
    sides.append(nodes[1].to_edge_array(nodes[2]))

    # area is half the size of the cross product
    if signed:
        return np.linalg.norm(np.cross(sides[0], sides[1]))/2.0

    return np.cross(sides[0], sides[1])/2.0

def edges_to_area_ratio_squared(nodes):
    """
    find the ratio of the square of the total edge length to the area: L^2/A
    """
    sides = []
    sides.append(nodes[0].to_edge_array(nodes[1]))
    sides.append(nodes[1].to_edge_array(nodes[2]))
    sides.append(nodes[2].to_edge_array(nodes[0]))

    # area is half the size of the cross product
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
        offsets.append(v3.Vector3([tmp_x, tmp_y, tmp_z]))

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

def dispersity(values):
    """
    find the variance of the values devided by the square of the average
    Args:
        values (numbers)
    Return:
        (numpy.float64): variance/average^2
    """
    return np.var(values)/(np.mean(values)**2)
