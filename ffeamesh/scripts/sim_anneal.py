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
import sys
import pathlib
import argparse

import ffeamesh.tetmeshtools.tetgenread as tr
import ffeamesh.tetmeshtools.ffeavolfilereader as fr
import ffeamesh.optimizemesh.costfunction as cf
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

    raise ValueError(f"probability {text} not in range [1.0, 0.0]")

def parse_weights(text0, text1):
    print(f"PW {text0} {text1}")
    raise ValueError("W")

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

def main():
    """run the script"""
    args = get_args()

    weights = None
    if args.weights is not None:
        print(args.weights)
        if len(args.weights) != 5:
            print(f"Wrong number of cost function weights {len(args.weights)}, should be five", file=sys.stderr)
            sys.exit(1)

        tmp = [x>=0.0 for x in args.weights]
        if not all(tmp):
            print(f"At least one cost function weight was negative, {args.weights}", file=sys.stderr)
            sys.exit(1)

        weights = cf.CostFeatures(inv_vol_dispersity = args.weights[0],
                                  shape_tets = args.weights[1],
                                  faces_shape = args.weights[2],
                                  total_shape = args.weights[3],
                                  isovalue_fit = args.weights[4])
    else:
        weights = cf.CostFeatures(inv_vol_dispersity = 1.0,
                                  shape_tets = 1.0,
                                  faces_shape = 1.0,
                                  total_shape = 1.0,
                                  isovalue_fit = 0.0)

    try:
        model = None
        if args.mesh_file.suffix == ".vol":
            model = fr.make_model_from_ffea(args.mesh_file)
        else:
            model = tr.make_model_from_tetgen(args.mesh_file)

        mutate = sa.MutateParams(args.mutate_prob, args.max_mutate, args.all_xyz_mutate)

        sa.simulated_anneal(args.cooling,
                            model,
                            weights,
                            mutate,
                            args.debug,
                            args.out_file_root)

    except (ValueError, IOError) as err:
        print(f"Error: {err}")
        return

if __name__ == "__main__":
    main()
