"""
Created on 19 April 2023

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this file except in compliance with the License. You may obtain a copy of the
License at http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

This work was funded by Joanna Leng's EPSRC funded RSE Fellowship (EP/R025819/1)

@copyright 2023
@author: j.h.pickering@leeds.ac.uk and j.leng@leeds.ac.uk
"""
# set up linting conditions
# pylint: disable = import-error
import pathlib
import argparse

import ffeamesh.optimizemesh.simanneal as sa

def handel_input_file(filename):
    """
    check input file and and convert to pathlib.Path
    Args:
        filename (str): input file
    Returns
        pathlib.Path
    """
    if filename.endswith(".vol"):
        return pathlib.Path(filename)

    tetgen_suffix = [".1.ele", ".1.node", ".1.face"]
    if filename.endswith(tetgen_suffix[0]):
        filename = filename.strip(tetgen_suffix[0])
    elif filename.endswith(tetgen_suffix[1]):
        filename = filename.strip(tetgen_suffix[1])
    elif filename.endswith(tetgen_suffix[2]):
        filename = filename.strip(tetgen_suffix[2])

    return pathlib.Path(filename)

def positive_int(text):
    """convert text to int and test >0"""
    value = int(text)
    if value > 0:
        return value

    raise ValueError(f"{text} is not >0")

def positive_float(text):
    """convert text to int and test >=0"""
    value = float(text)
    if value >= 0.0:
        return value

    raise ValueError(f"{text} is not >=0")

def probability(text):
    """convert text to a probabilty"""
    value = float(text)
    if 1.0 >= value >= 0.0:
        return value

    raise ValueError(f"probability {text} not in range [0.0, 0.1]")

def get_args():
    """
    get command line arguments
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("-d",
                        "--debug",
                        type=sa.DebugLevel,
                        default='none',
                        help="set debug level")

    parser.add_argument("-c",
                        "--cooling",
                        type=sa.CoolingFunction,
                        required=True,
                        choices=list(sa.CoolingFunction),
                        help="simulated anneal selection cooling function")

    parser.add_argument("-m",
                       "--mesh_file",
                       type=handel_input_file,
                       required=True,
                       help="name root of tetgen files, or ffea .vol file")

    parser.add_argument("-o",
                       "--out_file_root",
                       type=str,
                       required=False,
                       help="name root of output tetgen files")

    parser.add_argument("-p",
                       "--mutate_prob",
                       type=probability,
                       default=0.1,
                       help="probability that a node will be mutated")

    parser.add_argument("-x",
                       "--max_mutate",
                       type=positive_float,
                       default=0.1,
                       help="maximum size of mutation on a single axis")

    parser.add_argument("-a",
                       "--all_xyz_mutate",
                       action='store_true',
                       help="if true mutation applied to all three axis, else rnd selection")

    parser.add_argument("-w",
                        "--weights",
                        type=float,
                        nargs='+',
                        help="enter weights")

    return parser.parse_args()

def get_preprocess_args():
    """
    get command line arguments
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("-d",
                        "--debug",
                        type=sa.DebugLevel,
                        default='none',
                        help="set debug level")

    parser.add_argument("-m",
                       "--mesh_file",
                       type=handel_input_file,
                       required=True,
                       help="name root of tetgen files, or ffea .vol file")

    parser.add_argument("-o",
                       "--out_file_root",
                       type=str,
                       required=False,
                       help="name root of output tetgen files")

    parser.add_argument("-p",
                       "--mutate_prob",
                       type=probability,
                       default=0.1,
                       help="probability that a node will be mutated")

    parser.add_argument("-x",
                       "--max_mutate",
                       type=positive_float,
                       default=0.1,
                       help="maximum size of mutation on a single axis")

    parser.add_argument("-a",
                       "--all_xyz_mutate",
                       action='store_true',
                       help="if true mutation applied to all three axis, else rnd selection")

    parser.add_argument("-w",
                        "--weights",
                        type=float,
                        nargs='+',
                        required=True,
                        help="enter weights")

    parser.add_argument("-s",
                       "--steps",
                       type=positive_int,
                       default=100,
                       help="number of steps of random walk")

    return parser.parse_args()
