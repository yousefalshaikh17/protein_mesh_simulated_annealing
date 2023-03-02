
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
import collections

Node = collections.namedtuple("Node", "index, x, y, z")

def get_data_start_end_indices(lines, keyword):
    """
    get the start and end indices for a block beginning with keyword
    Args:
        lines (str): source text, newlines stripped
        keywork (str): text marker for start of block
    Returns:
        tuple(int, int): the start and end indices of data
    """
    key_line = lines.index(keyword)
    number_entries = int(lines[key_line+1])
    start = key_line+2

    return (start, start+number_entries)

def get_points(lines):
    """
    extract the points from the file contents
    Args:
        lines ([str]): text file lines with endlines removed
    Returns
        [[float, float, float]]: the points
    """
    start, end = get_data_start_end_indices(lines, "points")
    points = []
    for index in range(start, end):
        parts = lines[index].split()
        points.append([float(x) for x in parts])

    return points

def get_surface(lines):
    """
    extract the triangles from the file contents
    Args:
        lines ([str]): text file lines with endlines removed
    Returns
        [[int, int, int]]: the points
    """
    start, end = get_data_start_end_indices(lines, "surfaceelements")
    triangles = []
    for index in range(start, end):
        parts = lines[index].split()
        triangles.append([int(parts[5]), int(parts[6]), int(parts[7])])

    return triangles

def get_volume(lines):
    """
    extract the tets from the file contents
    Args:
        lines ([str]): text file lines with endlines removed
    Returns
        [[int, int, int, int]]: the points
    """
    start, end = get_data_start_end_indices(lines, "volumeelements")
    tets = []
    for index in range(start, end):
        parts = lines[index].split()
        tets.append([int(parts[2]), int(parts[3]), int(parts[4]), int(parts[5])])

    return tets

def read_file(source_file):
    """
    read ffea source file
    Args:
        source_file (pathlib.Path): the file
    Returns
        [[float, float, float]]: vertices
        [[int, int, int]]: surface triangles
        [[int, int, int, int]]: volume tets
    """
    with source_file.open(mode='r') as in_file:
        lines = in_file.readlines()

    lines = [x.strip() for x in lines]

    points = get_points(lines)
    surface = get_surface(lines)
    volume = get_volume(lines)

    return points, surface, volume
