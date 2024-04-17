Tetrahedral Meshing Tools
=========================

These Python tools were originally developed to support tetrahedral meshes for the Fluctuating Finite Element Analysis [FFEA](https://bitbucket.org/FFEA/ffea/downloads/) package, the authors gratefully acknowledge the help and support of the FFEA team.

The purpose is to convert Medical Research Council (MRC) files to tetrahedral mesh files in Tetgen format, view and analyse the results. FFEA output meshes in their own (.vol) format are also viewable.

The main tools comprise:

    * an converter from MRC files, which are used to store cryogenic electron microscope (cryo-em) data, to tetrahedral mesh files, in Tetgen or VTK formats.

    * a visualization tool that both allows the user to interactively explore the tetrahedral mesh and provides properties of the individual tetrahedra.

    * tools for extracting data from MRC files and thresholding the files.

This work was funded by Joanna Leng's EPSRC funded RSE Fellowship (EP/R025819/1)

Project start in summer 2021.

## Version

This is version 1.0

## Copyright and License

The algorithm and software in this project were developed by Joanna Leng, Jonathan Pickering and J Rogers together with the FFEA team at the University of Leeds. The main funding for this was Joanna Leng's Research Software Engineering Fellowship (EP/R025819/1).

Licensed under GNU General Public License v3.0 as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. You may not use this file except in compliance with the License. You may obtain a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.html>.

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

## Developed With

The project was developed using Python 3.9 and Anaconda, Inc. on Windows 10 systems. PyQt5 was and pyopengl were used in the viewer.

## Table of Contents

**[Quick Start](#quick-start)**<br>
**[Installation](#installation)**<br>
***[Python Environment](#python-environment)**<br>
***[Scripts](#scripts)**<br>
**[Usage](#usage)**<br>
**[Prerequisites](#prerequisites)**<br>
**[Installation](#istallation)**<br>
**[Unit Tests](#unit-tests)**<br>
**[Worked Example](#example)**<br>
**[Coding Converntions](#coding)**<br>

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

`git clone https://github.com/jonathanHuwP/tet_mesh_tools`

or, if you alreay have cloned

`git pull`

If you have an old conda environment remove it.

`conda env remove --name tetmesht`

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

## Scripts

The following scrips are available, all will provide instructions if run from the command line with either '-h' or '--help' options:

### MRC Data

   * **mrc_header_info.py** - print the header from an MRC file

   * **mrc_image_stats.py** - reads an MRC file and prints out its image intensity stats

   * **mrc_voxel_size.py** - print an MRC file's voxel size

   * **pdb_to_mrc.py** - simulate a cryo-em MRC file from a protein data base (PDB) using Van der Walls spheres.

### Filter

   * **mrc_crop.py** - crops 3D MRC image file data

   * **mrc_threshold.py** - reads an MRC file and set all voxels less than the threshold to zero, then outputs to a new file

   * **mrc_coarsen.py** - Coarsens MRC files to a user-defined resolution and outputs them as a new MRC file.

### Meshing

   * **mrc_to_tets.py** - process MRC files and produces a regular tetrahedral volumetric meshes using the "marching tetrahedra" algorithm. Decomposition of the image voxels into 5 or 6 tetrahedra each is available. Output can be in the Tetgen .ele, .face, and .node file format and/or the .vtk format.

### Viewing

   * **tgv.py** - run a viewer capable of rendering a three dimensional images of a mesh. Both Tetgen format files (.node, .ele, .face) and FFEA format files (.vol) can be read. Properties of the tetrahedra and overall mesh are also displayed. View modes are: the whole mesh; the outer surface as a wire frame; and any individual tetrahedron.

### Make Example Mesh Files

   * **make_tet_mesh_examples.py** - Make a pair of voxels decomposed into tetrahedra and output in VTK format. 5 and 6 tetrahedra decompositions are available.

## Usage

### mrc_to_tets

five_tets & six_tets input MRC files and processes them to produce volumetric mesh
 files. VTK files can also be produced for mesh analysis, using vtk based tools such as ParaView.

To get started quickly, the usage message can be viewed with the following command:

      python mrc_to_tets.py -h

In order to run five or six_tets there are three required flags:

   * -h, --help            show this help message and exit

   * `-i`, `--input` input file
   * `-o`, `--output` output file name root, (no suffix)

   * `-v`, `--vtk`     output files in vtk format.

   * `-f`, `--ftetg`    output files in Tetgen format.

   * `-t` THRESHOLD, `--threshold` lower filter for voxels, default zero

   * `-w`, `--overwrite` overwrite preexisting output files of same name

   * `-V`, `--verbose`         write verbose output

   * `-p`, `--progress`        print progress during operation (may slow package)

   * `-6`, `--use_six_tets`    decompose into six tets, default 5

   * `-m {1,2,3,4}`, `--low_vertices {1,2,3,4}`
               number vertices above isovalue for tet to be included in mesh (default: 2)

   * `-n VOX_COUNTS X Y Z`, `--vox_counts VOX_COUNTS X Y Z`
               voxel size for tets, if not used same as image, enter x y & z in Angstroms

The `-n` option allows the coarseness of the mesh to be specified, so `-n 5 10 15` will produce a mesh with five voxels on the x axis, ten on the y and fifteen on the z axis.

The `-m` option specifies the maximum number of vertices that can be below the isovalue on a tetrahedron that is included in the mesh, specified by the PruneLevel enumeration. For example, if prune level is set to 2 then tetrahedra with 0, 1 or 2 vertices below the isovalue are includes in the mesh, while tets with 3 or 4 vertices below the isovalue are culled.

||0|1|2|3|4|
|----|---|---|---|---|---|
|**PruneLevel.1**|In|In|**Out**|**Out**|**Out**|
|**PruneLevel.2**|In|In|In|**Out**|**Out**|
|**PruneLevel.3**|In|In|In|In|**Out**|
|**PruneLevel.4**|In|In|In|In|In|

**Table 1.** The rows specify the possible prune levels the user can choose. For each level the in the mesh, out of the mesh choice is specified for the five possible numbers of vertices a tetrahedron can have below the isovlue, from none to all four.

In the output directory defined by `-o` or `--output`, the following files will be found after running tet_from_pix successfully, where [outputname] is the user-defined name of the output files:

   * [outputname].1.ele
   * [outputname].1.face
   * [outputname].1.node
   * [outputname].vtk

## Prerequisites

   * [Python (>= 3.8)](https://www.python.org/).
     Required for running the tool scripts.

   * [NumPy (<1.23.0)](https://numpy.org/).
     Required Python library.

   * [SciPy](https://scipy.org/).
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

## Worked Example

Because publicly available PDB files are more common than cryo-em MRC files, in this example a MRC is simulated from a PDB. The example uses the small plant protein crambin (3nir.pdb) a available from [https://www.rcsb.org/structure/3NIR](https://www.rcsb.org/structure/3NIR).

1. Make a simulated MRC file by running `pdb_to_mrc -i <path>\3nir.pdb -o <path>\3nir.mrc -n 15 15 15 -w soft`

2. To see the information in the file's header run `mrc_header_info -i <path>\3nir.mrc`

3. To see the statistics of the data in the image run `mrc_image_stats -i <path>\3nir.mrc`

4. To convert to a mesh in Tetgen format using five tetrahedra per voxel, run: `mrc_to_tets -i <path>\3nir.mrc -o 3nir -v -f -t 1.5 -w -V -p -m2 -n 10 10 10`

5. View the mesh run `tgv -i <path>\3nir`

## Coding Converntions

The main convention we observe it that all code should have a pylint score above 7.5.
