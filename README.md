Protein Mesh Optimization Tool Using Simulated Annealing
=========================

This Python tool allows for a protein mesh to be optimized for use in simulation. Currently, the tool only optimizes the faces of the tetrahedron mesh and the fitness of each node.

The main tools comprise:

    * Read & Write functionality for tetrahedal meshes.

    * Ability to compute densities of mesh nodes, faces, and elements using an MRC file.

    * Perform Simulated Annealing on the mesh to improve its quality.

This work was developed for the Advanced Computer Science Masters of Science Degree Project at the University of Leeds. It makes use of the tet_mesh_tools repository which could be found here: <https://github.com/jonathanHuwP/tet_mesh_tools>

Project start in summer 2024.

## Version

This is version 1.0 which is what was provided for the project.

## Copyright and License

Licenses of the forked repositories apply.

The algorithm and software in the tet_mesh_tools repository were developed by Joanna Leng, Jonathan Pickering and J Rogers together with the FFEA team at the University of Leeds. The main funding for this was Joanna Leng's Research Software Engineering Fellowship (EP/R025819/1).

Licensed under GNU General Public License v3.0 as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. You may not use this file except in compliance with the License. You may obtain a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.html>.

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

## Developed With

The project was developed using Python 3.9 and Anaconda, Inc. on Windows 11 systems. tet_mesh_tools played a big role in the development of the software.

## Quick Start

This software runs in a Python 3.9 environment within Anaconda, so you will need to install and open an Anaconda shell. Once that is open, move to the top directory of tet_mesh_tools (the directory with the file README.md in it) and type the following the FIRST time you run the tet_mesh_tools software. Not all the instructions are required for later runs:

`conda env create -f environment_tetmesht.yml`

Next, activate the cpt Anaconda environment using the following command:

`conda activate tetmesht`

Now install the tools into the environment:

`pip install --editable . -v`

Execution of the visualizer:

`tgv`

Execution of a command line script is in this format when you download the mvc file called filename.

`mrc_header_info -i filename`

And to get help on how to use a script:

`mrc_header_info --help`

If you want to use them to write your own python scripts you can now import them into a script.

`import tetmeshtools.<import name>....`

### Qt User Interface

The user interface of the viewer is written in Qt, using PyQt5. To change the main window itself edit the `resources/tetgenviewermain.ui` file. Editing is best done using QtDesigner, which can be run from the command line by typing `designer`. If callback function names are changed then manual editing is required.

To produce the Ui_tetgenviewermain.py file you have to use the PyQt pyuic5 program.

'pyuic5 .\resources\tetgenviewermain.ui -o .\tetmeshtools\app_tgv\gui\Ui_tetgenviewermain.py'

## Installation

This software was written to be run in a conda environment.

### Python Environment:

We advise you to get your conda installation up to date before you do this, but this is not obligatory with the command:

`conda update --all`

To see what conda environments you have, run the command.

`conda env list`

### Installing the Tools

Get a copy of the latest version from GitHub by running.

`git clone https://github.com/yousefalshaikh17/protein_mesh_simulated_annealing`

or, if you alreay have cloned

`git pull`

If you have an old conda environment (or the tet_mesh_tools environment) remove it.

`conda env remove --name tetmesht`

CD into the repository.

`cd protein_mesh_simulated_annealing`

Create a new environment.

`conda env create -f .\environment_tetmesht.yml`

Activate the environment.

`conda activate tetmesht`

Install the package to the environment using pip.

`pip install --editable .`

This will produce a directory called 'tetmeshtools.egg-info', do not move, rename, or delete this directory.

Make the documentation run doxygen

`doxygen`

To stop using the enviroment:

`conda deactivate`

To remove the environment:

`conda env remove --name tetmesht`

### Building the Documentation

The documentation is build using doxygen,and can be built by running the command 'doxygen' in the command tool window. This will create a directory 'doc/html' holding the documentation web pages, the root page being 'index.html'. To modify doxygen's behaviour edit the file 'Doxyfile' either manually or using doxywizard. For instance, to add Latex output change the flag on line 1806 from 'NO' to 'YES'.

## Scripts

The following scrips are available, all will provide instructions if run from the command line with either '-h' or '--help' options:

### Mesh Management

   * **copy_mesh.py** - Makes an unmodified copy of the input mesh file.

### MRC-Mesh Tools

   * **compute_densities.py** - Prints density values of mesh nodes, faces, and tetrahedron given MRC Image and mesh.

### Optimization

   * **simulated_annealing.py** - Optimizes a mesh given an MRC Image and mesh.

## Usage

### simulated_annealing

This script performs simulated annealing on a mesh given the MRC image, the tetrahedral mesh, and other crucial parameters.

To get started quickly, the usage message can be viewed with the following command:

      python simulated_annealing.py -h

In order to run five or six_tets there are three required flags:

   * -h, --help            show this help message and exit

   * `-m`, `--mrcfile` input MRC file (required)

   * `-i`, `--tetmeshfile` input tetmesh file (required)

   * `-o`, `--output` output file name root, (no suffix) (required)

   * `-A`, `--density_weight` sum of density weights, default 1

   * `-B`, `--surface_area_weights` weights for variance of triangle surface area, default 1

   * `--mutation_probability` probability of mutating a node (required)

   * `--mutation_multiplier` multiplier controlling how large mutation changes can be (required)

   * `-t`, `--temperature`, `--initial_temperature` initial temperature for the annealing process (required)

   * `--temperature_decrement` decrement of temperature, default 1

   * `--iterations_per_temperature` number of iterations per temperature change (required)

   * `-s`, `--seed` seed to control randomness

   * `--target_density_sum` target density sum during annealing (required)

   * `-k` multiplier against temperature for checking probability of keeping worse fitness, default 1

   * `-p`, `--progress` print optimization progress (may slow down optimization)

In the output directory defined by `-o` or `--output`, the optimized mesh can be found after running the script successfully, where [outputname] is the user-defined name of the output files:

   * [outputname].1.ele
   * [outputname].1.face
   * [outputname].1.node

## Prerequisites

   * [Python (>= 3.8)](https://www.python.org/).
     Required for running the tool scripts.

   * [NumPy (<1.23.0)](https://numpy.org/).
     Required Python library.

   * [SciPy](https://scipy.org/).
     Required Python library.
     
   * [pandas](https://pandas.pydata.org/).
     Required Python library.

   * [VTK for Python](https://pypi.org/project/vtk/).
     Required Python library.

   * [mrcfile for Python](https://pypi.org/project/mrcfile/).
     Required Python library.

   * [doxygen](https://pypi.org/project/doxypypy/)
     Required for building documentation

   * [pylint](https://pypi.org/project/pylint/)
     Required for development

## Unit Tests

The directory `tests` contain unit tests for developers working on tet_mesh_tools. Inside is a script called `unit_tests.py` which executes the tests, this can be run by typing.

`python tetmeshtools/tests/unit_tests.py`

### Unit Tests Added

   * Test to ensure that exported mesh is identical to the original mesh.
   * Test to ensure that a node is mutated within the limits of the mutation.
   * Test to ensure that the probability of acceptance works as intended for the Metropolis Algorithm.
