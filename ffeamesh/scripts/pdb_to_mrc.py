
#!/usr/bin/env python
"""
 mrc_coarsen.py

 A script that produces lower resolution mrc files.

 ----------------------------

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
# set up linting
# pylint: disable = import-error

import argparse
import pathlib
import typing
import itertools
import sys
import getpass
import datetime
import time
import multiprocessing as multi
import numpy as np

import Bio.PDB

from ffeamesh.mrc_utility import (CellSize, CellAngles, CellProps, write_mrcfile)

## map of atomic radii (Wikipedia)
_atomic_rad ={'H': 1.2,
              'C': 1.7,
              'N': 1.6,
              'O': 1.5,
              'P': 1.8,
              'S': 1.8}

class AtomBall(typing.NamedTuple):
    """storage for an atom as radius and location"""
    radius: float
    x     : float
    y     : float
    z     : float

    def __str__(self):
        """string reper for user"""
        return f"<AtomBall, {self.radius}, ({self.x}, {self.y}, {self.z})>"

class SpaceBounds():
    """spatial bounds of collection of points"""

    def __init__(self, x_coord, y_coord, z_coord):
        """
        set up object
        Args:
            x_coord (float): x coordinate of first point
            y_coord (float): y coordinate of first point
            z_coord (float): z coordinate of first point
        """
        ## low x bound
        self.low_x = x_coord
        ## high x bound
        self.high_x = x_coord

        ## low y bound
        self.low_y = y_coord
        ## high y bound
        self.high_y = y_coord

        ## low z bound
        self.low_z = z_coord
        ## high z bound
        self.high_z = z_coord

    def add_point(self, x_coord, y_coord, z_coord):
        """
        add a point to the collection
        Args:
            x_coord (float): x coordinate of first point
            y_coord (float): y coordinate of first point
            z_coord (float): z coordinate of first point
        """
        if x_coord < self.low_x:
            self.low_x = x_coord
        elif x_coord > self.high_x:
            self.high_x = x_coord

        if y_coord < self.low_y:
            self.low_y = y_coord
        elif y_coord > self.high_y:
            self.high_y = y_coord

        if z_coord < self.low_z:
            self.low_z = z_coord
        elif z_coord > self.high_z:
            self.high_z = z_coord

    def pad(self, padding=5.0):
        """
        add padding to limits
        Args:
            padding (float): the padding
        """
        self.low_x  -= padding
        self.high_x += padding

        self.low_y  -= padding
        self.high_y += padding

        self.low_z  -= padding
        self.high_z += padding

    def get_cella(self):
        """
        get mrc cell dimensions
        Returns:
            mrc_utility.CellSize
        """
        return CellSize(self.high_x-self.low_x,
                        self.high_y-self.low_y,
                        self.high_z-self.low_z)

    def get_origin(self):
        """
        get mrc origin
        Returns:
            [float, float, float]: x, y, z
        """
        return [self.low_x, self.low_y, self.low_z]

    def __str__(self):
        """
        user description of object
        Returns
            string
        """
        str_x = f"x({self.low_x}, {self.high_x})"
        str_y = f"y({self.low_y}, {self.high_y})"
        str_z = f"z({self.low_z}, {self.high_z})"
        return f"<SpaceBounds: {str_x}, {str_y}, {str_z}"

class VoxelSizes(typing.NamedTuple):
    """storage for an voxel size and half sizes"""

    ## size size x
    size_x: np.float32

    ## half size size x
    half_size_x: np.float32

    ## size size y
    size_y: np.float32

    ## half size size y
    half_size_y: np.float32

    ## size size z
    size_z: np.float32

    ## half size size z
    half_size_z: np.float32

class VoxelModel():
    """
    mapping from indices to space
    """

    def __init__(self, bounds, num_x, num_y, num_z):
        """
        setup object
            bounds (SpaceBounds): spacial extent
            num_x (int): number voxels on x
            num_y (int): number voxels on y
            num_z (int): number voxels on z
        """
        ## pointer to the bounds object
        self.bounds = bounds

        ## number of voxels on x axis
        self.num_x = num_x

        ## number of voxels on y axis
        self.num_y = num_y

        ## number of voxels on z axis
        self.num_z = num_z

        # voxel sizes
        size_x = bounds.get_cella().x/num_x
        size_y = bounds.get_cella().y/num_y
        size_z = bounds.get_cella().z/num_z

        ## storage for the size sizes
        self.sizes = VoxelSizes(size_x, size_x/2, size_y, size_y/2, size_z, size_z/2)

    def to_cartesian(self, index_x, index_y, index_z):
        """
        convert voxel indices to cartesian coordinates of centre
        Args:
            index_x (int):
            index_y (int):
            index_z (int):
        Returns
            (float, float, float): cartesian coordinatess
        """
        coord_x = self.bounds.low_x + (index_x * self.sizes.size_x) + self.sizes.half_size_x
        coord_y = self.bounds.low_y + (index_y * self.sizes.size_y) + self.sizes.half_size_y
        coord_z = self.bounds.low_z + (index_z * self.sizes.size_z) + self.sizes.half_size_z

        return (coord_x, coord_y, coord_z)

    def shape(self):
        """
        get shape of array
        Reurns:
            (int, int, int): couns in (z, y, x)
        """
        return (self.num_z, self.num_y, self.num_x)

    def number_of_voxels(self):
        """
        total number of voxels
        Returns:
            int
        """
        return self.num_x*self.num_y*self.num_z

    def __str__(self):
        """
        to text rep
        """
        counts = f"({self.num_x}, {self.num_y}, {self.num_z})"
        sizes = f"({self.sizes.size_x:.3f}, {self.sizes.size_y:.3f}, {self.sizes.size_z:.3f})"

        return f"<VoxelModel: counts {counts}, size {sizes}>"

def make_bounds(atoms):
    """
    get a bound object from a list of atoms
    Args:
        atoms [AtomBall]
    returns:
        SpaceBounds
    """
    bounds = SpaceBounds(atoms[0].x, atoms[0].y, atoms[0].z)

    for atom in atoms[1:]:
        bounds.add_point(atom.x, atom.y, atom.z)

    bounds.pad()
    return bounds

def atom_radius(element):
    """
    get radius for element
    Args:
        element (str)
    Return
        float: atomic radius Angstoms
    """
    if element in _atomic_rad:
        return _atomic_rad[element]

    return 1.5

def read_pdb(file_path):
    """
    read contents of pdb file
    Args:
        file_path (pathlib.Path):
    Returns
        Bio.SeqRecord
    Raises:
        Exception if file cannot be read
    """
    pdbparser = Bio.PDB.PDBParser(QUIET=True)   # suppress PDBConstructionWarning
    struct = pdbparser.get_structure(None, file_path)

    atoms = []
    for atom in struct.get_atoms():
        ball = AtomBall(atom_radius(atom.element),
                        atom.coord[0],
                        atom.coord[1],
                        atom.coord[2])
        atoms.append(ball)

    return atoms

def density_at(atoms, model, index_x, index_y, index_z):
    """
    make the mrc data field
    Args:
        atoms [AtomBall]:
        model VoxelModel:
        index_x int: array index of voxel
        index_y int: array index of voxel
        index_z int: array index of voxel
    Returns:
        float density at ctr of voxel
    """
    voxel_ctr = model.to_cartesian(index_x, index_y, index_z)
    density = 0.0

    for atom in atoms:
        distance = np.linalg.norm([voxel_ctr[0] - atom.x,
                                   voxel_ctr[1] - atom.y,
                                   voxel_ctr[2] - atom.z])

        twice_rad = 2*atom.radius
        if distance < atom.radius:
            density += 1
        elif distance < twice_rad:
            density += (twice_rad-distance)/atom.radius

    return density

def output_details(input_file, atoms, bounds):
    """
    print out details of file
    Args:
        input_file (pathlib.Path): source file
        atoms [AtomBall]: list of atoms
        bounds VoxelBounds: bounds for atoms
    """
    print(f"Source: {input_file}", file=sys.stdout)
    print(f"\tNumber of atoms: {len(atoms)}", file=sys.stdout)

    del_x = bounds.high_x - bounds.low_x
    del_y = bounds.high_y - bounds.low_y
    del_z = bounds.high_z - bounds.low_z

    print(f"\tSize of cell (x, y, z): ({del_x:.3f}, {del_y:.3f}, {del_z:.3f})")

def make_mrc_data(atoms, model, proc_num, print_progress=True):
    """
    make the mrc data array
    Args:
        atoms [AtomBall]:
        model VoxelModel:
        proc_num int: the process number
        print_progress bool: if true print every 1000 iterations
    """
    data = np.zeros(model.shape(), dtype=np.float32)
    total = model.number_of_voxels()
    count = itertools.count(1)
    for index_z in range(model.num_z):
        for index_y in range(model.num_y):
            for index_x in range(model.num_x):
                data[index_z, index_y, index_x] = density_at(atoms,
                                                             model,
                                                             index_x,
                                                             index_y,
                                                             index_z)
                done = next(count)
                if print_progress and (done%1000 == 0):
                    print(f"\tProcess {proc_num}: has completed {done} out of {total} voxels", file=sys.stdout)

    if print_progress:
        print(f"\tProcess {proc_num} finished: {total} voxels completed", file=sys.stdout)

    return data

def run_multiprocess(atoms, model):
    """
    run the density calculation in parallel
    Args:
        atoms [AtomBall]
        model [VoxelModel]
    Returns:
        numpy.ndarray: the voxel densities
    """
    num_processors = multi.cpu_count()
    if num_processors>2:
        num_processors -= 1

    chunk = int(len(atoms)/num_processors)
    pool = multi.Pool(num_processors)
    processes = []

    t0 = time.time()
    # distribute the atoms array across the processors
    for count in range(num_processors):
        start = count*chunk
        if count == num_processors-1:
            # final case must go to end of array to allow for rounding error in size of chunk
            #data.append(make_mrc_data(atoms[start:], model))
            proc_args = (atoms[start:], model, count)
        else:
            #data.append(make_mrc_data(atoms[start:(start+chunk)], model))
            proc_args = (atoms[start:(start+chunk)], model, count)

        processes.append(pool.apply_async(make_mrc_data, proc_args))

    data = [p.get() for p in processes]
    # add the densities in the arrays
    total = data[0]
    for tmp in data[1:]:
        total = np.add(total, tmp)

    print(f"Data completed elapsed time {time.time()-t0:.2f}s", file=sys.stdout)

    return total

def write_out_file(data, bounds, in_file, out_file):
    """
    write out the MRC file
    Args:
        data [AtomBall]
        bounds SpaceBounds
        in_file pathlib.Path
        out_file pathlib.Path
    """
    label = f"From {str(in_file)}"
    label += f" by {getpass.getuser()} on "
    label += datetime.datetime.now().strftime("%d-%b-%Y (%H:%M:%S.%f)")

    angles = CellAngles(90.0, 90.0, 90.0)
    write_mrcfile(data,
                  CellProps(bounds.get_cella(), angles),
                  out_file,
                  label,
                  out_file,
                  origin=bounds.get_origin())

    print(f"Data written to {out_file}", file=sys.stdout)

def pdb_path_exists(file_name):
    """
    convert str to pathlib.Path and check existance
    Args
        file_name (str)
    Returns
        pathlib.Path
    Raises
        ValueError if non-existing
    """
    path = pathlib.Path(file_name)

    if not path.exists():
        raise ValueError(f"File {file_name} does not exist!")

    if path.suffix.lower() != '.pdb':
        raise ValueError(f"File {file_name} does not end with .pdb")

    return path

def get_args():
    """
    get the command line arguments
    """
    description = """make a simulated MRC format electron density map from a
PDB format list of atom locations, densities are calculated
at the centre point of each voxel bases on VDW radii of the atoms"""

    parser = argparse.ArgumentParser(description)

    parser.add_argument('-i',
                        '--input',
                        type=pdb_path_exists,
                        required=True,
                        help='source pdb file (must end .pdb or .PDB)')

    parser.add_argument('-o',
                        '--output',
                        type=pathlib.Path,
                        default=pathlib.Path('pdb_to_mrc_OUT.mrc'),
                        help='source pdb file')

    parser.add_argument("-n",
                        '--num_voxels',
                        type=int,
                        nargs='+',
                        help='number of voxels on x y z axis')

    parser.add_argument("-w",
                        "--overwrite",
                        action="store_true",
                        help="if output file exists overwrite")

    return parser.parse_args()

def run_pdb_to_mrc():
    """
    get pdb file from user an run pdb to mrc on it
    """
    args = get_args()

    # validate number of voxels
    if args.num_voxels is not None:
        if len(args.num_voxels) != 3:
            sys.exit("you must provide three counts to the num_voxels (run with -h for help)")
        if not all(x > 0 for x in args.num_voxels):
            sys.exit("voxel counts must be greater than zero")

    # read file, process atoms and print results
    atoms = read_pdb(args.input)
    bounds = make_bounds(atoms)
    output_details(args.input, atoms, bounds)

    if args.num_voxels is None:
        return

    model = VoxelModel(bounds, args.num_voxels[0], args.num_voxels[1], args.num_voxels[2])
    data = run_multiprocess(atoms, model)
    write_out_file(data, bounds, args.input, args.output)

if __name__ == "__main__":
    run_pdb_to_mrc()
