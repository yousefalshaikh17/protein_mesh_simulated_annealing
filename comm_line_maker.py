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
"""
import argparse
import pathlib
import csv
import collections

ParameterSet = collections.namedtuple("ParameterSet", [
                                      "name",
                                      "density",
                                      "shear_visc",
                                      "bulk_visc",
                                      "shear_mod",
                                      "bulk_mod",
                                      "stokes_rad"])

def get_args():
    """
    get the command line arguments
        Returns
            (argparse.namespace)
    """
    parser = argparse.ArgumentParser("""make command line for running voltoffea""")

    parser.add_argument("-i",
                        "--input",
                        type=pathlib.Path,
                        required=True,
                        help="input csv file")

    parser.add_argument("-r",
                        "--root_name",
                        type=str,
                        required=True,
                        help="root name of ele, face & node files")

    parser.add_argument("-n",
                        "--set_number",
                        type=int,
                        required=True,
                        help="set number of parameters to be used")

    return parser.parse_args()

def main():
    """
    run the program
    """
    args = get_args()
    parameters = read_params(args.input)

    make_command(args.root_name, parameters[args.set_number])

def read_params(csv_file):
    """read the parameters"""

    if not csv_file.exists():
        raise ValueError(f"""File "{csv_file}" does not exist""")

    with csv_file.open('r') as file_in:
        reader = csv.reader(file_in)

        # remove headers
        next(reader)

        parameters = []
        for row in reader:
            parameters.append(ParameterSet(row[0].strip(),
                                           row[1].strip(),
                                           row[2].strip(),
                                           row[3].strip(),
                                           row[4].strip(),
                                           row[5].strip(),
                                           row[6].strip()))

        return parameters

def make_command(root_name, params):
    """
    write out the new command
    Args
        root_name (str): root name of input files
        params (ParameterSet): the parameter values
    """
    command = []
    command.append("ffeatools voltoffea")
    command.append(f"""--mesh {root_name+".1.ele"} {root_name+".1.face"} {root_name+".1.node"}""")
    command.append(f"--density {params.density}")
    command.append(f"--shear-visc {params.shear_visc}")
    command.append(f"--bulk-visc {params.bulk_visc}")
    command.append(f"--shear-mod {params.shear_mod}")
    command.append(f"--bulk-mod {params.bulk_mod}")
    command.append(f"--stokes-radius {params.stokes_rad}")
    command.append("--make-script")

    print(' '.join(command))

if __name__ == "__main__":
    main()
