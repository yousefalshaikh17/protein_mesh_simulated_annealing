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
resolution = 20
scale_factor=20/(mrc.voxel_size.tolist()[0])
b = scipy.ndimage.zoom(a, scale_factor, order=3)
with mrcfile.new('data/whole_0p5.mrc', overwrite=True) as mrc_2:
    mrc_2.set_data(b)
    mrc_2.header.origin = mrc.header.origin
    new_vox = mrc.voxel_size.tolist()
    new_list = []
    for i in new_vox:
        new_list.append((1/scale_factor)*i)
    mrc_2.voxel_size = tuple(new_list)
    
# print(np.shape(b))

new_mrc= mrcfile.open('data/whole_0p5.mrc', mode='r+')
print(new_mrc.header)

mrcfile.validate('data/whole_0p5.mrc')