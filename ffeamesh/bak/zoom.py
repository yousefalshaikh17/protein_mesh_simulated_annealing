#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 11:50:39 2021

@author: mollygravett
"""

import mrcfile
import numpy as np
import scipy
from scipy import ndimage, misc
mrc=mrcfile.open('data/map_equilibrium_mystructrigor_15A_0p00202.mrc', mode='r+')
a = mrc.data
nx = mrc.header.nx
ny = mrc.header.ny
nz = mrc.header.nz

resolution = 20
scale_factor=(mrc.voxel_size.tolist()[0])/resolution


print(scale_factor)
b = scipy.ndimage.zoom(a, scale_factor)

with mrcfile.new('data/whole_20_correct_error_real.mrc', overwrite=True) as mrc_2:
    mrc_2.set_data(b)
    mrc_2.header.origin = mrc.header.origin
    old_vox = mrc.voxel_size.tolist()
    new_list =[] 
    nx_2 = mrc_2.header.nx
    ny_2 = mrc_2.header.ny
    nz_2 = mrc_2.header.nz
    scale_x = nx_2/nx
    scale_y = ny_2/ny
    scale_z = nz_2/nz
    new_list.append((1/scale_x)*old_vox[0]*(nx_2/(nx_2-1)))
    new_list.append((1/scale_y)*old_vox[1]*(ny_2/(ny_2-1)))
    new_list.append((1/scale_z)*old_vox[2]*(nz_2/(nz_2-1)))
    mrc_2.voxel_size = tuple(new_list)

#not sure if doing different scale for each dimension necessary. maybe fine to just use scale_x for all?

