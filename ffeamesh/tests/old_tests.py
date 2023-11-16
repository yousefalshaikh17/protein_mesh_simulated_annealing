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
        unit_tests.py
        Author: Jarvellis Rogers - University of Leeds
        Email: J.F.Rogers1@leeds.ac.uk
"""

import getopt
import mrcfile
import numpy as np
import subprocess
import sys

def bash_cmd(cmd):
    subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, executable="/bin/bash")

def usage():
    helpMessage = """   Coding:   Jarvellis Rogers (J.F.Rogers1@leeds.ac.uk)

Unit tests for ffea-meshing.

General Options:
  -h [ --help ]             Print usage message.
  -t [ --test ] (=1,2,3)    Run specific unit tests.

Unit Tests:
  1) Tests running tet_from_pix without coarsening the input MRC file.
  2) Tests running mrc_zoom.
  3) Tests running tet_from_pix with coarsening of the input MRC file.
"""

    print(helpMessage)
    sys.exit()

try:
    options, remainder = getopt.getopt(sys.argv[1:], "ht:", ["help", "test="])
except getopt.GetoptError as err:
    print("ERROR: " + str(err) + "\n")

test1 = False
test2 = False
test3 = False

for opt, arg in options:
    if opt in ("-h", "--help"):
        usage()
    elif opt in ("-t", "--test"):
        try:
            tests = list(map(int, arg.split(",")))
        except ValueError:
            msg = "ERROR: Value for tests must be 1, 2, and/or 3."
            sys.exit(msg)
        for i in tests:
            if i == 1:
                test1 = True
            elif i == 2:
                test2 = True
            elif i == 3:
                test3 = True
            else:
                msg = "ERROR: Value for tests must be 1, 2, and/or 3."
                sys.exit(msg)

# Turns tests on if none were specified in the command line
if test1 == False and test2 == False and test3 == False:
    test1 = True
    test2 = True
    test3 = True

if test1:
    print("Running Test 1: tet_from_pix without coarsening")
    bash_cmd("python ../tet_from_pix.py --input=input.mrc --output=output/test1_result --threshold=0.00202")

    # Open tet files for diffing
    with open("test1_output.1.ele") as ele1:
        ele1Data = ele1.read()
    with open("output/test1_result.1.ele") as ele2:
        ele2Data = ele2.read()
    with open("test1_output.1.face") as face1:
        face1Data = face1.read()
    with open("output/test1_result.1.face") as face2:
        face2Data = face2.read()
    with open("test1_output.1.node") as node1:
        node1Data = node1.read()
    with open("output/test1_result.1.node") as node2:
        node2Data = node2.read()
    with open("test1_output.vtk") as vtk1:
        vtk1Data = vtk1.read()
    with open("output/test1_result.vtk") as vtk2:
        vtk2Data = vtk2.read()

    # Diffs original output files and test result files
    testPass = True
    if ele1Data == ele2Data:
        print("\t.ele files match")
    else:
        print("\tWARNING: .ele files do not match")
        testPass = False
    if face1Data == face2Data:
        print("\t.face files match")
    else:
        print("\tWARNING: .face files do not match")
        testPass = False
    if node1Data == node2Data:
        print("\t.node files match")
    else:
        print("\tWARNING: .node files do not match")
        testPass = False
    if vtk1Data == vtk2Data:
        print("\t.vtk files match")
    else:
        print("\tWARNING: .vtk files do not match")
        testPass = False

    if testPass:
        print("Test 1: Passed")
    else:
        print("Test 1: Failed")

if test2:
    print("Running Test 2: mrc_zoom")
    bash_cmd("python ../mrc_zoom.py --input=input.mrc --output=output/test2_result --resolution=20")

    # Open MRC files for comparrison
    mrc1 = mrcfile.open("test2_output.mrc", mode='r')
    mrc1Data = mrc1.data
    mrc2 = mrcfile.open("output/test2_result.mrc", mode='r')
    mrc2Data = mrc2.data

    # Compares MRC files
    testPass = True
    if np.array_equal(mrc1Data, mrc2Data):
        print("\t.mrc files match")
    else:
        print("\tWARNING: .mrc files do not match")
        testPass = False

    if testPass:
        print("Test 2: Passed")
    else:
        print("Test 2: Failed")

if test3:
    print("Running Test 3: tet_from_pix with coarsening")
    bash_cmd("python ../tet_from_pix.py --input=input.mrc --output=output/test3_result --threshold=0.00202 --resolution=20")

    # Open tet files for diffing
    with open("test3_output.1.ele") as ele1:
        ele1Data = ele1.read()
    with open("output/test3_result.1.ele") as ele2:
        ele2Data = ele2.read()
    with open("test3_output.1.face") as face1:
        face1Data = face1.read()
    with open("output/test3_result.1.face") as face2:
        face2Data = face2.read()
    with open("test3_output.1.node") as node1:
        node1Data = node1.read()
    with open("output/test3_result.1.node") as node2:
        node2Data = node2.read()
    with open("test3_output.vtk") as vtk1:
        vtk1Data = vtk1.read()
    with open("output/test3_result.vtk") as vtk2:
        vtk2Data = vtk2.read()

    # Open MRC files for comparrison
    mrc1 = mrcfile.open("test3_output.mrc", mode='r')
    mrc1Data = mrc1.data
    mrc2 = mrcfile.open("output/test3_result.mrc", mode='r')
    mrc2Data = mrc2.data

    testPass = True

    # Diffs original output files and test result files
    if ele1Data == ele2Data:
        print("\t.ele files match")
    else:
        print("\tWARNING: .ele files do not match")
        testPass = False
    if face1Data == face2Data:
        print("\t.face files match")
    else:
        print("\tWARNING: .face files do not match")
        testPass = False
    if node1Data == node2Data:
        print("\t.node files match")
    else:
        print("\tWARNING: .node files do not match")
        testPass = False
    if vtk1Data == vtk2Data:
        print("\t.vtk files match")
    else:
        print("\tWARNING: .vtk files do not match")
        testPass = False

    # Compares MRC files
    if np.array_equal(mrc1Data, mrc2Data):
        print("\t.mrc files match")
    else:
        print("\tWARNING: .mrc files do not match")
        testPass = False

    if testPass:
        print("Test 3: Passed")
    else:
        print("Test 3: Failed")
