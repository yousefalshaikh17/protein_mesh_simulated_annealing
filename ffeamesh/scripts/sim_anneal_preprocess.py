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
import mrcfile

import ffeamesh.tetmeshtools.tetgenread as tr
import ffeamesh.tetmeshtools.ffeavolfilereader as fr
import ffeamesh.optimizemesh.costfunction as cf
import ffeamesh.optimizemesh.simanneal as sa
import ffeamesh.scripts.simannealcomline as cl
import ffeamesh.mrclininterpolate as mi

def make_weights(weights_vec):
    """
    converts weights vector to a CostFratures object
    Args:
        weights_vec (list(float))
    Returns:
        CostFeatures
    """
    length = len(weights_vec)
    if length != 6:
        print(f"Wrong number of cost function weights: {length}, should be six",
                file=sys.stderr)
        sys.exit(1)

    tmp = [x>=0.0 for x in weights_vec]
    if not all(tmp):
        print(f"At least one cost function weight was negative, {weights_vec}",
                file=sys.stderr)
        sys.exit(1)

    return  cf.CostFeatures(inv_vol_dispersity = weights_vec[0],
                            shape_tets = weights_vec[1],
                            faces_shape = weights_vec[2],
                            total_shape = weights_vec[3],
                            isovalue_fit = weights_vec[4],
                            dist_to_image2= weights_vec[5])

def main():
    """run the script"""
    args = cl.get_preprocess_args()

    weights = make_weights(args.weights)

    if not args.mrc_file.exists():
        print(f"Error MRC file {args.mrc_file} does not exist.")
        return

    with mrcfile.open(args.mrc_file, mode='r+') as mrc:
        image = mi.MRCImage(mrc)

        try:
            model = None
            if args.mesh_file.suffix == ".vol":
                model = fr.make_model_from_ffea(args.mesh_file)
            else:
                model = tr.make_model_from_tetgen(args.mesh_file)

            mutate = sa.MutateParams(args.mutate_prob, args.max_mutate, args.all_xyz_mutate)

            sa.random_walk(model,
                        weights,
                        args.isovalue,
                        image,
                        mutate,
                        args.steps,
                        args.debug)

        except (ValueError, IOError) as err:
            print(f"Error: {err}")

if __name__ == "__main__":
    main()
