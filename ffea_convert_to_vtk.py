# -*- coding: utf-8 -*-
#
#  This file is part of the FFEA simulation package
#
#  Copyright (c) by the Theory and Development FFEA teams,
#  as they appear in the README.md file.
#
#  FFEA is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  FFEA is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
#
#  To help us fund FFEA development, we humbly ask that you cite
#  the research papers on the package.
#

"""
        FFEA_convert_to_VTK.py
        Author: Jarvellis Rogers, University of Leeds
        Email: j.f.rogers1@leeds.ac.uk
        Requirements: Python >= 3.8
"""

import getopt
import sys
import re
from dataclasses import dataclass
from os import path
from typing import List, Optional, Tuple
from scipy.spatial.transform import Rotation as rot
import numpy as np

def path_check(fpath, i_o):
    if not path.exists(fpath):
        ERROR_MESSAGE = "ERROR: " + i_o + " file path does not exist. Please try again."
        sys.exit(ERROR_MESSAGE)

try:
    options, remainder = getopt.getopt(sys.argv[1:], "i:o:", ["input=", "output="])
except getopt.GetoptError as err:
    print("ERROR: " + str(err) + "\n")

for opt, arg in options:
    if opt in ("-i", "--input"):
        ffeaFilePath = arg
        path_check(ffeaFilePath, "Input")
    elif opt in ("-o", "--output"):
        destFilePath = arg
        path_check(destFilePath, "Output")

try:
    ffeaFilePath
except NameError:
    ERROR_MESSAGE = "ERROR: FFEA file path not defined."
    sys.exit(ERROR_MESSAGE)
try:
    destFilePath
except NameError:
    ERROR_MESSAGE = "ERROR: Log files destination path not defined."
    sys.exit(ERROR_MESSAGE)

inputFilePath = path.dirname(ffeaFilePath)

@dataclass
class Object:
    nodes: str = ""
    centroid: Tuple[int,int,int] = (0,0,0)
    rotation: Tuple[float,float,float,float,float,float,float,float,float] = (1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0) # Rotation matrix of no rotation
    scale: float = float()

def centroid_parse(x):
    x = x.replace("(","")
    x = x.replace(")","")
    x = tuple(map(float, x.split(',')))
    return x

def rotation_parse(x):
    x = x.replace("(","")
    x = x.replace(")","")
    x = list(map(float, x.split(',')))
    if len(x) == 3:
        x = rot.from_euler("xyz",x,degrees=True)
        x = x.as_matrix()
        y = []
        for i in x:
            for j in i:
                y.append(j)
        return y
    elif len(x) == 9:
        return x
    else:
        print("WARNING: One of the rotation values in your FFEA script appears to be neither a Euler angle or rotation matrix. This run will probably be unsucessful.")
        return x

# Regex to match against
blobOpen = re.compile(r"<\s*blob\s*>")
blobClose = re.compile(r"<\/\s*blob\s*>")
rodOpen = re.compile(r"<\s*rod\s*>")
rodClose = re.compile(r"<\/\s*rod\s*>")
nodeTag = re.compile(r"<\s*nodes\s*=\s*([^>]+)>")
inputTag = re.compile(r"<\s*input\s*=\s*([^>]+)>")
centroidTag = re.compile(r"<\s*centroid\s*=\s*([^>]+)>")
centroidPosTag = re.compile(r"<\s*centroid_pos\s*=\s*([^>]+)>")
rotationTag = re.compile(r"<\s*rotation\s*=\s*([^>]+)>")
scaleTag = re.compile(r"<\s*scale\s*=\s*([^>]+)>")
numNodes = re.compile(r"num_nodes\s*[^>\n]+")
nodeLine = re.compile(r"\s*[+-]?\d+\.?\d*\s*[+-]?\d+\.?\d*\s*[+-]?\d+\.?\d*")
numElements = re.compile(r"num_elements\s*,[^>\n]+")
elementsPrevLine = re.compile(r"FRAME\s*\d+\s*ROD\s*\d+")

OPEN_SEARCH = True # Checks if we're searching for new blob or rod group or if one has been found
BLOB_OPEN_TAG = False
ROD_OPEN_TAG = False
FIRST_TAG_FOUND = False
CURR_OBJECT: Optional[Object] = None
objects: List[Object] = []

# Parses FFEA scripts and appends nodes and centroids to a Blob dataclass list
with open(ffeaFilePath,"r") as f:
    for ln in f:
        ln = ln.strip()
        if OPEN_SEARCH:
            if blobOpen.match(ln):
                OPEN_SEARCH = False
                BLOB_OPEN_TAG = True
                continue
            elif rodOpen.match(ln):
                OPEN_SEARCH = False
                ROD_OPEN_TAG = True
                continue

        elif BLOB_OPEN_TAG:
            if (tag := nodeTag.match(ln)):
                if FIRST_TAG_FOUND:
                    assert CURR_OBJECT
                    CURR_OBJECT.nodes = tag.group(1)
                else:
                    assert not CURR_OBJECT
                    CURR_OBJECT = Object(nodes=tag.group(1))
                    FIRST_TAG_FOUND = True
                continue
            if (tag := centroidTag.match(ln)):
                if FIRST_TAG_FOUND:
                    assert CURR_OBJECT
                    x = tag.group(1)
                    x = centroid_parse(x)
                    CURR_OBJECT.centroid = x
                else:
                    assert not CURR_OBJECT
                    x = tag.group(1)
                    x = centroid_parse(x)
                    CURR_OBJECT = Object(centroid=x)
                    FIRST_TAG_FOUND = True
                continue
            if (tag := rotationTag.match(ln)):
                if FIRST_TAG_FOUND:
                    assert CURR_OBJECT
                    x = tag.group(1)
                    x = rotation_parse(x)
                    CURR_OBJECT.rotation = x
                    FIRST_TAG_FOUND = True
                else:
                    assert not CURR_OBJECT
                    x = tag.group(1)
                    x = rotation_parse(x)
                    CURR_OBJECT = Object(rotation=x)
                continue
            if (tag := scaleTag.match(ln)):
                if FIRST_TAG_FOUND:
                    assert CURR_OBJECT
                    CURR_OBJECT.scale = float(tag.group(1))
                else:
                    assert not CURR_OBJECT
                    CURR_OBJECT = Object(scale=float(tag.group(1)))
                    FIRST_TAG_FOUND = True
                continue
            if blobClose.match(ln):
                objects.append(CURR_OBJECT)
                CURR_OBJECT = None
                OPEN_SEARCH = True
                FIRST_TAG_FOUND = False
                BLOB_OPEN_TAG = False
                continue

        elif ROD_OPEN_TAG:
            if (tag := inputTag.match(ln)):
                if FIRST_TAG_FOUND:
                    assert CURR_OBJECT
                    CURR_OBJECT.nodes = tag.group(1)
                else:
                    assert not CURR_OBJECT
                    CURR_OBJECT = Object(nodes=tag.group(1))
                    FIRST_TAG_FOUND = True
                continue
            if (tag := centroidPosTag.match(ln)):
                if FIRST_TAG_FOUND:
                    assert CURR_OBJECT
                    x = tag.group(1)
                    x = centroid_parse(x)
                    x = tuple([float(float(i) * 1.7e-10) for i in x]) # Centroids need to be converted to rod units in order to apply them - https://ffea.readthedocs.io/en/stable/namespacemesoDimensions.html#a2c51b401ba77d7948a6e153599451679
                    CURR_OBJECT.centroid = x
                else:
                    assert not CURR_OBJECT
                    x = tag.group(1)
                    x = centroid_parse(x)
                    x = tuple([float(float(i) * 1.7e-10) for i in x])
                    CURR_OBJECT = Object(centroid=x)
                    FIRST_TAG_FOUND = True
                continue
            if (tag := rotationTag.match(ln)):
                if FIRST_TAG_FOUND:
                    assert CURR_OBJECT
                    x = tag.group(1)
                    x = rotation_parse(x)
                    CURR_OBJECT.rotation = x
                    FIRST_TAG_FOUND = True
                else:
                    assert not CURR_OBJECT
                    x = tag.group(1)
                    x = rotation_parse(x)
                    CURR_OBJECT = Object(rotation=x)
                continue
            if (tag := scaleTag.match(ln)):
                if FIRST_TAG_FOUND:
                    assert CURR_OBJECT
                    CURR_OBJECT.scale = float(tag.group(1))
                else:
                    assert not CURR_OBJECT
                    CURR_OBJECT = Object(scale=float(tag.group(1)))
                    FIRST_TAG_FOUND = True
                continue
            if rodClose.match(ln):
                objects.append(CURR_OBJECT)
                CURR_OBJECT = None
                OPEN_SEARCH = True
                FIRST_TAG_FOUND = False
                ROD_OPEN_TAG = False
                continue
# print("DEBUG - Objects:")
# print(objects) #DEBUG

nodeCoords = []
ffeaDir = ffeaFilePath[:ffeaFilePath.rindex("/")+1]
totalNodes = [] # For counting the number of nodes to add to the VTK file

# Parses object coordinates and applies any translations
for i in objects:
    CURR_TOTAL_NODES = 0
    nodeDir = ffeaDir + i.nodes
    currCoords = []
    ELEMENTS_PREV_FOUND = False # Checking to see if elem coords are found in .rod files

    # Parses .node and .rod files and extracts all coordinates
    with open(nodeDir,"r") as f:
        if nodeDir.endswith(".node"):
            for ln in f:
                if numNodes.match(ln):
                    CURR_TOTAL_NODES = re.compile(r'(\d+)$').search(ln).group(1)
                if nodeLine.match(ln):
                    lineCoords = []
                    x = [float(j) for j in ln.split()]
                    lineCoords.append(x)
                    currCoords.append(lineCoords)
        elif nodeDir.endswith(".rod"):
            for ln in f:
                if numElements.match(ln):
                    CURR_TOTAL_NODES = re.compile(r'(\d+)$').search(ln).group(1)
                if elementsPrevLine.match(ln):
                    ELEMENTS_PREV_FOUND = True
                    continue
                if ELEMENTS_PREV_FOUND:
                    lineCoords = []
                    x = [float(j) for j in ln.split(",")]
                    lineCoords.append(x)
                    groupCoords = np.array(lineCoords).reshape(-1,3)
                    for k in groupCoords:
                        currCoords.append(k)
                    ELEMENTS_PREV_FOUND = False

    if int(CURR_TOTAL_NODES) != len(currCoords):
        print("WARNING: The number of nodes/elements found does not match the stated amount in " + i.nodes + ". Something may have gone wrong in your file.")
        print("Number of nodes stated by " + i.nodes + ": " + str(CURR_TOTAL_NODES))
        print("Number of nodes found: " + str(len(currCoords)))

    x = np.array(i.centroid)
    y = np.array(currCoords)

    # Applies rotation
    y = y.flatten()
    y = np.array(y).reshape(-1,3)
    r = np.array(i.rotation)
    r = np.array(r).reshape(-1,3)
    r = rot.from_matrix(r)
    y = r.apply(y)

    # Blobs' default centroids are the average of node coords, needs to be subtracted to get it to (0,0,0) like rods
    # Many thanks to Molly Gravett for figuring this headache out
    if nodeDir.endswith(".node"):
        yAvg = np.average(y, axis=0)
        y = y - yAvg

    # Applies scaling
    if i.scale != float(0):
        if nodeDir.endswith(".node") and i.scale != 1e-10:
            scaleFactor = i.scale / 1e-10 # Gets ratio between recorded scale and the default 1e-10
            x = x * scaleFactor
            y = y * scaleFactor
        elif nodeDir.endswith(".rod") and i.scale != 1:
            x = x * i.scale
            y = y * i.scale

    # Applies centroid
    z = x + y

    # Converts rod units to angstroms for consistency with blobs
    if nodeDir.endswith(".rod"):
        z = z * 1e10

    nodeCoords.append(z)
    totalNodes.append(int(CURR_TOTAL_NODES))

# Sets up the VTK file and header information
totalNodes = sum(totalNodes)
ffeaFilename = ffeaFilePath[ffeaFilePath.rindex("/")+1:]
ffeaFilename = ffeaFilename[:ffeaFilename.index(".")]
vtkFilename = ffeaFilename + ".vtk"
vtkFilepath = destFilePath + "/" + vtkFilename

# Creates a VTK file and appends the node coordinates
with open(vtkFilepath, 'w') as f:
    f.write("# vtk DataFile Version 5.1\n")
    f.write("vtk file generated by FFEA_convert_to_VTK.py\n")
    f.write("ASCII\n")
    f.write("DATASET UNSTRUCTURED_GRID\n")
    f.write("POINTS " + str(totalNodes) + " float\n")

    for i in nodeCoords:
        for j in i:
            j = str(j)
            j = j.replace("[","")
            j = j.replace("]","")
            j = j.strip()
            f.write(j + "\n")

#TODO: Coupling
# in x to y, y is what moves position
