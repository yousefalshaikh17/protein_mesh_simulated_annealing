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

    Authors: Joanna Leng, Jonathan Pickering - University of Leeds
    Emails: J.Leng@leeds.ac.uk, J.H.Pickering@leeds.ac.uk
"""
import numpy as np
from itertools import product, count
import ffeamesh.voxels2tets_utility as v2t
import sys

class Grid():
    """
    make and store a 3D grid
    """

    def __init__(self, counts, start, end, image_counts):
        """
        set up object
        Args:
            counts ([int, int, int]): number of steps on x, y & z axis
            start ([float, float, float]): minimum values of x, y & z
            end ([float, float, float]): maximum values of x, y & z
            image_counts ([int, int, int]): number of voxels on x, y & z in image
        """
        ## number of steps on x, y & z axis
        self._num_steps = counts

        ## mimimum coordinates of image cube
        self._start = [np.float64(x) for x in start]
        ## maximum coordinates of image cube
        self._end = [np.float64(x) for x in end]

        ## offsets to allow for half voxel bound on linear interpolation
        self._offsets = []
        for index in range(3):
            self._offsets.append((self._end[index]-self._start[index])/(2*image_counts[index]))

        # move start & end in by half an image voxel
        for index in range(3):
            self._start[index] = self._start[index]+self._offsets[index]
            self._end[index] = self._end[index]-self._offsets[index]

        ## flat list of all 3D vertices
        self._vertices = None

        ## storage for the desnsities at each vertex
        self._densities = []

        ## storage for the connectivit of the tets
        self._connectivities = []

    def build_grid(self, image):
        """
        build grid and calculate densities
        Args:
            image (MRCImage): the density interpolater
        """
        # make points on each axis
        lin_x = np.linspace(self._start[0],
                            self._end[0],
                            self._num_steps[0]+1)

        lin_y = np.linspace(self._start[1],
                            self._end[1],
                            self._num_steps[1]+1)

        lin_z = np.linspace(self._start[2],
                            self._end[2],
                            self._num_steps[2]+1)

        # outer product of axis points to make 3D
        tmp = list(product(lin_z, lin_y))
        tmp = list(product(tmp, lin_x))

        # flat list of all 3D vertices
        self._vertices = [[x[1], x[0][1], x[0][0]] for x in tmp]

        self._calculate_densities(image)

    def _calculate_densities(self, mrc):
        """
        find & store image density for each vertex
        """
        for point in self._vertices:
            density, distance = mrc.density_or_distance_at(point[0], point[1], point[2])
            if distance != 0.0:
                raise ValueError(f"mesh point outside image {distance}")
            self._densities.append(density)

    def build_five_tets(self):
        """
        make the tet's connectivity
        """
        # iterate the voxels and make tet connectivities
        for voxel_z in range(self._num_steps[2]):
            for voxel_y in range(self._num_steps[1]):
                for voxel_x in range(self._num_steps[0]):

                    if v2t.is_odd(voxel_x, voxel_y, voxel_z):
                        indices = v2t.odd_cube_tet_indices()
                    else:
                        indices = v2t.even_cube_tet_indices()

                    verts_indices = self.get_voxel_indices(voxel_x, voxel_y, voxel_z)

                    # make the tets
                    for tet in indices:
                        tet_indices = []
                        for index in tet:
                            # need to add the indices in the original array return array of indices
                            tet_indices.append(verts_indices[index])
                        self._connectivities.append(tet_indices)

    def build_six_tets(self):
        """
        make the tet's connectivity
        TODO complete
        """
        # iterate the voxels and make tet connectivities
        for voxel_z in range(self._num_steps[2]):
            for voxel_y in range(self._num_steps[1]):
                for voxel_x in range(self._num_steps[0]):

                    # if v2t.is_odd(voxel_x, voxel_y, voxel_z):
                    #     indices = v2t.odd_cube_tet_indices()
                    # else:
                    #     indices = v2t.even_cube_tet_indices()
                    indices = v2t.cube_6_tet_indices()

                    verts_indices = self.get_voxel_indices(voxel_x, voxel_y, voxel_z)

                    # # make the tets
                    for tet in indices:
                        tet_indices = []
                        for index in tet:
                            # need to add the indices in the original array return array of indices
                            tet_indices.append(verts_indices[index])
                        self._connectivities.append(tet_indices)

    def crop_mesh_to_isovalue(self, isovalue, prune_level, progress):
        """
        Remove tets outside or largly outside isovalue
        Args:
            isovalue (float): limit value
            level (PruneLevel): the number of vertices below the isovalue that causes deleation
            progress (bool): if true print progress reports
        """
        completed = count()
        total_tets = len(self._connectivities)
        new_tets = []

        for tet in self._connectivities:
            count_outside = 0
            for index in tet:
                if self._densities[index] < isovalue:
                    count_outside += 1

            if count_outside <= prune_level.value:
                new_tets.append(tet)

            tmp = next(completed)
            if progress and tmp%10000 == 0:
                print(f"Processed {tmp} out of {total_tets} tets", file=sys.stdout)

        if progress:
            print(f"Crop of tets completed: {len(new_tets)} surviving tets")

        self._connectivities = new_tets

    def remove_surplas_vertices(self):
        """
        remove the vertices not used in the connectivity and relabel the connectivity
        """
        used = set()

        for tet in self._connectivities:
            for index in tet:
                used.add(index)

        saved_verts = []
        map_to_new_indices = {}
        new_index = count()
        for index in used:
            saved_verts.append(self._vertices[index])
            map_to_new_indices[index] = next(new_index)

        new_connectivities = [[None, None, None, None] for _ in self._connectivities]

        for tet_index, tet in enumerate(new_connectivities):
            for vert_index in range(4):
                old_index = self._connectivities[tet_index][vert_index]
                tet[vert_index] = map_to_new_indices[old_index]

        self._connectivities = new_connectivities
        self._vertices = saved_verts

    def get_total_num_voxels(self):
        """
        get the total number of voxels
        Return:
            int
        """
        return self._num_steps[0]*self._num_steps[1]*self._num_steps[2]

    def get_num_voxels_x(self):
        """
        get number of points on x axis
        Return:
            int: the number of points
        """
        return self._num_steps[0]

    def get_num_voxels_y(self):
        """
        get number of points on y axis
        Return:
            int: the number of points
        """
        return self._num_steps[1]

    def get_num_voxels_z(self):
        """
        get number of points on z axis
        Return:
            int: the number of points
        """
        return self._num_steps[2]

    def get_num_verts_x(self):
        """
        get number of points on x axis
        Return:
            int: the number of points
        """
        return self._num_steps[0]+1

    def get_num_verts_y(self):
        """
        get number of points on y axis
        Return:
            int: the number of points
        """
        return self._num_steps[1]+1

    def get_num_verts_z(self):
        """
        get number of points on z axis
        Return:
            int: the number of points
        """
        return self._num_steps[2]+1

    def get_vertex_index(self, x_index, y_index, z_index):
        """
        get the indices into the vertices array for the given 3D indices
        Args:
            x_index (int): the index of point on x axis
            y_index (int): the index of point on y axis
            z_index (int): the index of point on z axis
        Return:
            [float, float, float]: the point
        """
        z_offset = z_index * (self._num_steps[1]+1) * (self._num_steps[0]+1)
        y_offset = y_index * (self._num_steps[0]+1)
        return z_offset + y_offset + x_index

    def get_vertex(self, x_index, y_index, z_index):
        """
        get point
        Args:
            x_index (int): the index of point on x axis
            y_index (int): the index of point on y axis
            z_index (int): the index of point on z axis
        Return:
            [float, float, float]: the point
        """
        return self._vertices[self.get_vertex_index(x_index, y_index, z_index)]

    def get_vertices(self):
        """
        get all points
        Return:
            [[float, float, float]..]: the points
        """
        return self._vertices

    def get_densities(self):
        """
        get all densities
        Return:
            [float]: the densities
        """
        return self._densities

    def get_connectivities(self):
        """
        get all densities
        Return:
            [float]: the densities
        """
        return self._connectivities

    def get_voxel_indices(self, x_index, y_index, z_index):
        """
        get the array indices of the vertices froming the corners of a voxel
        Args:
            x_index (int): the x index of the voxel
            y_index (int): the y index of the voxel
            z_index (int): the z index of the voxel
                  Vertex Indices:
         7+----------+6
         /|         /|
        / |        / |
      4+----------+5 |
       |  |       |  |         Axes:
       | 3+-------|--+2        z  y
       | /        | /          | /
       |/         |/           |/
      0+----------+1           +----x
        Returns:
            [int x *]: indices of vertices froming the corners
        """
        indices = []
        indices.append(self.get_vertex_index(x_index, y_index, z_index))
        indices.append(self.get_vertex_index(x_index+1, y_index, z_index))
        indices.append(self.get_vertex_index(x_index+1, y_index+1, z_index))
        indices.append(self.get_vertex_index(x_index, y_index+1, z_index))

        indices.append(self.get_vertex_index(x_index, y_index, z_index+1))
        indices.append(self.get_vertex_index(x_index+1, y_index, z_index+1))
        indices.append(self.get_vertex_index(x_index+1, y_index+1, z_index+1))
        indices.append(self.get_vertex_index(x_index, y_index+1, z_index+1))

        return indices

    def get_voxel(self, x_index, y_index, z_index):
        """
        get the corners of a voxel
        Args:
            x_index (int): the x index of the voxel
            y_index (int): the y index of the voxel
            z_index (int): the z index of the voxel
        Returns:
            [[float, float, float]x8]: coordinates of corners
        """
        corners = []
        for index in self.get_voxel_indices(x_index, y_index, z_index):
            corners.append(self._vertices[index])

        return corners

    def get_start(self):
        """
        get mimimum coordinates of image cube
        Returns:
            [float, float, float]
        """
        return self._start

    def get_end(self):
        """
        get mimimum coordinates of image cube
        Returns:
            [float, float, float]
        """
        return self._end

    def __len__(self):
        """
        dunder method to allow objects to have len applied to them
        Return
            int: the lenght of the points array
        """
        return len(self._vertices)

    def __repr__(self):
        """
        dunder method for to string
        Returns:
            string
        """
        return f"Grid({self._num_steps}, {self._start}, {self._end})"
