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

import ffeamesh.tetmeshtools.tetgenread as tr
import ffeamesh.tetmeshtools.ffeavolfilereader as fr
import ffeamesh.optimizemesh.costfunction as cf
import ffeamesh.optimizemesh.simanneal as sa
import ffeamesh.scripts.simannealcomline as cl

def main():
    """run the script"""
    args = cl.get_args()

    weights = None
    if args.weights is not None:
        length = len(args.weights)
        if length != 5:
            print(f"Wrong number of cost function weights: {length}, should be five",
                  file=sys.stderr)
            sys.exit(1)

        tmp = [x>=0.0 for x in args.weights]
        if not all(tmp):
            print(f"At least one cost function weight was negative, {args.weights}",
                  file=sys.stderr)
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
