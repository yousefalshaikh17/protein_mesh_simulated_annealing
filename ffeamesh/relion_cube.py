#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 11:50:39 2021

@author: mollygravett
"""
#making into a box of equal dimensions for relion_image_handler

import mrcfile
import numpy as np

mrc=mrcfile.open('data/map_equilibrium_mystructrigor_15A_0p00202.mrc', mode='r+')
a = mrc.data

#making into cube
nxnynz = list(mrc.header.tolist()[0:3])
nz = nxnynz[2]
ny = nxnynz[1]
nx = nxnynz[0]
sorted_dimensions = sorted(nxnynz)
new_dimensions = (sorted_dimensions[-1],sorted_dimensions[-1],sorted_dimensions[-1])
new_dims_zeros = np.zeros(new_dimensions, dtype='float32')
for z in range(0, nz):
    for y in range(0, ny):
        for x in range(0, nx):
            new_dims_zeros[z, y, x] = a[z,y,x]

with mrcfile.new('data/box.mrc', overwrite=True) as mrc_box:
    mrc_box.set_data(new_dims_zeros)
    mrc_box.header.origin = mrc.header.origin
    mrc_box.voxel_size = mrc.voxel_size
