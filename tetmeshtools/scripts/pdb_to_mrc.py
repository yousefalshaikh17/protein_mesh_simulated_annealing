#!/usr/bin/env python
"""
 A script that produces lower resolution mrc files.

 ----------------------------

Licensed under the GNU General Public License, Version 3.0 (the "License"); you
may not use this file except in compliance with the License. You may obtain a
copy of the License at <https://www.gnu.org/licenses/gpl-3.0.html>.

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

This work was funded by Joanna Leng's EPSRC funded RSE Fellowship (EP/R025819/1)

@copyright 2023
@author: j.h.pickering@leeds.ac.uk and j.leng@leeds.ac.uk
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

from tetmeshtools.mrc_utility import (CellSize, CellAngles, CellProps, write_mrcfile)

## map of atomic radii (Wikipedia)
_atomic_rad ={'H': 1.2,
              'C': 1.7,
              'N': 1.6,
              'O': 1.5,
              'P': 1.8,
              'S': 1.8}

class Progress():
    """record progress in conversion"""

    def __init__(self, num_voxels, num_atoms):
        """
        set up object
        Args:
            num_voxels int:
            num_atoms int:
        """
        ## number of atom voxel calculations performed
        self._done = 0

        ## total numbe of atom voxel calculation
        self._total = num_voxels * num_atoms

    def __call__(self, done):
        """
        add a number of atom voxel calculations to running total
        Args:
            done int: numbe of atom voxel calculations performed
        Returns:
            float: current percentage of total
        """
        self._done += done
        return round(100.0*(self._done/self._total), 3)

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

def density_at(atoms, model, index_x, index_y, index_z, soft=None):
    """
    make the mrc data field
    Args:
        atoms [AtomBall]:
        model VoxelModel:
        index_x int: array index of voxel
        index_y int: array index of voxel
        index_z int: array index of voxel
        soft float: if not none exponential decay constant for soft atoms
    Returns:
        float density at ctr of voxel
    """
    voxel_ctr = model.to_cartesian(index_x, index_y, index_z)
    density = 0.0

    for atom in atoms:
        distance = np.linalg.norm([voxel_ctr[0] - atom.x,
                                   voxel_ctr[1] - atom.y,
                                   voxel_ctr[2] - atom.z])
        if distance < atom.radius:
            density += 1
        elif soft is not None:
            twice_rad = 2*atom.radius
            if distance < twice_rad:
                density += np.exp(soft*(distance-atom.radius))

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

def make_mrc_data(atoms, model, soft, proc_label, pipe):
    """
    make the mrc data array
    Args:
        atoms [AtomBall]:
        model VoxelModel:
        soft float: if not none exponential decay constant for soft atoms
        proc_label int: the process number
        pipe multiprocessing.Pipe: communication to parent
    """
    data = np.zeros(model.shape(), dtype=np.float32)
    count = itertools.count(1)

    for index_z in range(model.num_z):
        for index_y in range(model.num_y):
            for index_x in range(model.num_x):
                data[index_z, index_y, index_x] = density_at(atoms,
                                                             model,
                                                             index_x,
                                                             index_y,
                                                             index_z,
                                                             soft)

                done = next(count)
                if done%100 == 0:
                    # emit the number of atom calculations done
                    pipe.send((100*len(atoms), proc_label))

    pipe.send(('Fertig', (done%100)*len(atoms), proc_label))
    pipe.close()

    return data

def process_message(message, connections, progress):
    """
    handle a message
    Args:
        message (int/str, int):
        connections Connections:
        progress Progress: convert increment of operations to percentage of total
    """
    if isinstance(message[0], int):
        print(f"Conversion {progress(message[0]):.3f}% complete.", end='\r')
    else:
        percent = progress(message[1])
        print(f"Process {message[2]} finished: conversion {percent:.3f}% complete.")
        connections[message[2]].close()
        del connections[message[2]]

def setup_processes(pool, processes, pipes, model, atoms, soft, num_processors):
    """
    distribute the atoms array across the processors
    Args:
        pool multiprocessing.pool.Pool
        processes []
        pipes []
        model VoxelModel
        atoms [AtomBall]
        soft float: if not none exponential decay constant for soft atoms
        num_processors omt
    """
    chunk = int(len(atoms)/num_processors)

    for count in range(num_processors):
        start = count*chunk
        parent_conn, child_conn = multi.Pipe(False)
        if count == num_processors-1:
            # final case must go to end of array to allow for rounding error in size of chunk
            proc_args = (atoms[start:], model, soft, count, child_conn)
        else:
            proc_args = (atoms[start:(start+chunk)], model, soft, count, child_conn)

        processes.append(pool.apply_async(make_mrc_data, proc_args))
        pipes[count] = parent_conn

def manage_processes(pipes, progress):
    """
    handel messages from the processes, and close pipes as needed
    Args
        pipes [multiprocessing.Pipe]
        progress Progerss
    """
    flag = True
    while flag:
        for label in list(pipes.keys()):
            try:
                message = pipes[label].recv()
                process_message(message, pipes, progress)
            except EOFError:
                pipes[label].close()
                del pipes[label]
                print(f"Warning, process {label} exited abnormally!", file=sys.stderr)

            flag = bool(pipes) # False if empty

def run_pool(pool, model, atoms, soft, num_processors):
    """
    run the pool of processes
    Args:
        pool multiprocessing.pool.Pool
        model VoxelModel
        atoms [AtomBall]
        soft float: if not none exponential decay constant for soft atoms
        num_processors int
    """
    processes = []
    pipes = {}
    progress = Progress(model.number_of_voxels(), len(atoms))

    setup_processes(pool, processes, pipes, model, atoms, soft, num_processors)
    manage_processes(pipes, progress)

    return [p.get() for p in processes]

def run_multiprocess(atoms, model, soft=None):
    """
    run the density calculation in parallel
    Args:
        atoms [AtomBall]
        model [VoxelModel]
        soft float: if not none exponential decay constant for soft atoms
    Returns:
        numpy.ndarray: the voxel densities
    """
    num_processors = multi.cpu_count()
    if num_processors>2:
        num_processors -= 1

    time_start = time.time()
    with multi.Pool(num_processors) as pool:
        data = run_pool(pool, model, atoms, soft, num_processors)

    # add the densities in the arrays
    total = data[0]
    for tmp in data[1:]:
        total = np.add(total, tmp)

    print(f"Data completed elapsed time {time.time()-time_start:.2f}s", file=sys.stdout)

    return total

def make_mrc_label(in_file):
    """
    make array of label strings (max 80 character)
    Args:
        in_file (pathlib.Path)
    Returns:
        [str]
    """
    label = []

    label.append("Simulated MRC")

    label.append(f"{in_file.name[:80]}")
    if len(in_file.name) > 80:
        message = 'Warning file name truncated to 80 characters'
        label.append(message)
        print(message, file=sys.stderr)

    label.append(f"User: {getpass.getuser()}")

    now = datetime.datetime.now()
    label.append(now.strftime("Date: %d-%b-%Y"))
    label.append(now.strftime("Time: %H:%M:%S ") + now.astimezone().tzname())

    return label

def write_out_file(data, bounds, in_file, out_file):
    """
    write out the MRC file
    Args:
        data [AtomBall]
        bounds SpaceBounds
        in_file pathlib.Path
        out_file pathlib.Path
    """
    label = make_mrc_label(in_file)

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

def negative_float(text):
    """
    convert string to float and ensure negative
    Args:
        text str:
    Returns
        float < 0.0
    Raises
        ValueError
    """
    number = float(text)

    if number >= 0.0:
        raise ValueError("Exp decay constant must be less than zero")

    return number

def get_args():
    """
    get the command line arguments
    """
    description = """make a simulated MRC format electron density map from a
PDB format list of atom locations, densities are calculated
at the centre point of each voxel bases on VDW radii of the atoms"""

    parser = argparse.ArgumentParser(description)
    subparser = parser.add_subparsers(required=False)

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

    msg="use exponential decay of electron density outside VDW radii run soft -h for details"
    soft_p = subparser.add_parser("soft", help=msg)

    soft_p.add_argument('-e',
                        '--exp_const',
                        type=negative_float,
                        default = -6.0,
                        help='a in exp(a*r) for density outside VDW radius (default -6)')

    return parser.parse_args()

def main():
    """
    get pdb file from user an run pdb to mrc on it
    """
    args = get_args()
    soft_atoms = None
    if hasattr(args, "exp_const"):
        soft_atoms = args.exp_const

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
    data = run_multiprocess(atoms, model, soft_atoms)
    write_out_file(data, bounds, args.input, args.output)

if __name__ == "__main__":
    main()
