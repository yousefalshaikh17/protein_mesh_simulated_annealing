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
import csv

import ffeamesh.optimizemesh.simanneal as sa

class SimAnnealParameters():
    def __init__(self):
        self.debug = None
        self.cooling = None
        self.mesh_file = None
        self.out_file_root = None
        self.mutate_prob = None
        self.max_mutate = None
        self.all_xyz_mutate = None
        self.weights = None
        self.isovalue = None
        self.mrc_file = None
        self.start_temp = None
        self.stop_temp = None

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

def parse_weights(text):
    """convert text to weights"""
    return [float(part) for part in text.split()]

def params_from_str(params_dict):
    """
    convet paramters to their type
    """
    params = SimAnnealParameters()

    for key in params_dict:
        if key == "debug":
            params.debug = pathlib.Path(params_dict[key])
        elif key == "cooling":
            params.cooling = sa.CoolingFunction(params_dict[key])
        elif key == "mesh_file":
            params.mesh_file = handel_input_file(params_dict[key])
        elif key == "out_file_root":
            params.out_file_root = params_dict[key]
        elif key == "mutate_prob":
            params.mutate_prob = probability(params_dict[key])
        elif key == "max_mutate":
            params.max_mutate = positive_float(params_dict[key])
        elif key == "all_xyz_mutate":
            params.all_xyz_mutate = bool(params_dict[key])
        elif key == "weights":
            params.weights = parse_weights(params_dict[key])
        elif key == "isovalue":
            params.isovalue = positive_float(params_dict[key])
        elif key == "mrc_file":
            params.mrc_file = pathlib.Path(params_dict[key])
        elif key == "start_temp":
            params.start_temp = positive_float(params_dict[key])
        elif key == "stop_temp":
            params.stop_temp = positive_float(params_dict[key])
        elif key == "cooling_rate":
            params.cooling_rate = positive_float(params_dict[key])

    return params

def read_parameters(file):
    """
    get a dict of parameters from the
    Args:
        file (pathlib.Path): the parmeters file
    Retruns
        (dict): holding the parameters
    """
    params = {}
    with file.open('r') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='#')
        for row in reader:
            if len(row)>1:
                params[row[0]] = row[1].strip()

    return params_from_str(params)

def get_args():
    """
    get command line arguments
    """
    parser = argparse.ArgumentParser("""optimizes a tetrahedral mesh using simulated annealing""")
    subparser = parser.add_subparsers(required=True)

    sp_raw = subparser.add_parser("raw",
                                  help="input all parameter from command line")

    sp_raw.add_argument("-d",
                        "--debug",
                        type=pathlib.Path,
                        required=False,
                        help="save run data to debug file")

    sp_raw.add_argument("-c",
                        "--cooling",
                        type=sa.CoolingFunction,
                        required=True,
                        help=f"simulated anneal selection cooling function {[el.value for el in sa.CoolingFunction]}")

    sp_raw.add_argument("-r",
                        "--cooling_rate",
                        type=positive_float,
                        required=True,
                        help=f"cooling function rate (alpha) value")

    sp_raw.add_argument("-m",
                       "--mesh_file",
                       type=handel_input_file,
                       required=True,
                       help="name root of tetgen files, or ffea .vol file")

    sp_raw.add_argument("-o",
                       "--out_file_root",
                       type=str,
                       required=False,
                       help="name root of output tetgen files")

    sp_raw.add_argument("-p",
                       "--mutate_prob",
                       type=probability,
                       default=0.1,
                       help="probability that a node will be mutated")

    sp_raw.add_argument("-x",
                       "--max_mutate",
                       type=positive_float,
                       default=0.1,
                       help="maximum size of mutation on a single axis")

    sp_raw.add_argument("-a",
                       "--all_xyz_mutate",
                       action='store_true',
                       help="if true mutation applied to all three axis, else rnd selection")

    sp_raw.add_argument("-w",
                        "--weights",
                        type=float,
                        nargs='+',
                        help="enter weights")

    sp_raw.add_argument("-v",
                        "--isovalue",
                        type=float,
                        required=True,
                        help="isovalue defining surface in mrc file.")

    sp_raw.add_argument("-f",
                        "--mrc_file",
                        type=pathlib.Path,
                        required=True,
                        help="mrc image file on which isosurface is defined")

    sp_raw.add_argument("-start",
                        "--start_temp",
                        type=positive_float,
                        required=True,
                        help="start temperature for cooling schedule")

    sp_raw.add_argument("-stop",
                        "--stop_temp",
                        type=positive_float,
                        required=True,
                        help="stop temperature for cooling schedule")

    sp_file = subparser.add_parser("file",
                                   help="read parameters from a file")
    sp_file.set_defaults(read_func=read_parameters)

    sp_file.add_argument("-f",
                         "--file",
                         type=pathlib.Path,
                         required=True,
                         help="file of input parameters")

    return parser.parse_args()

def get_preprocess_args():
    """
    get command line arguments
    """
    parser = argparse.ArgumentParser()

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

    parser.add_argument("-v",
                        "--isovalue",
                        type=float,
                        required=True,
                        help="isovalue defining surface in mrc file.")

    parser.add_argument("-f",
                        "--mrc_file",
                        type=pathlib.Path,
                        required=True,
                        help="mrc image file on which isosurface is defined")

    return parser.parse_args()
