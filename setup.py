#!/usr/bin/env python
"""
setup.py

pip install the scripts

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
from setuptools import setup

DISTUTILS_DEBUG=1

setup(
    name='ffeamesh',
    version='1.0.0',
    author="The FFEA Team",
    description ='Meshing tools for the Fluctuating Finite Element Analysis tool (FFEA).',
    url='http://ffea.bitbucket.com',
    license='GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007',
    packages=['ffeamesh'],
    install_requires=[
        'argparse',
        'numpy',
        'vtk',
        'mrcfile'
    ],
    entry_points={
        'console_scripts': [
            'make_test_mrc_file = ffeamesh.scripts.make_test_mrcfile:main',
            'make_tet_mesh_examples = ffeamesh.scripts.make_tet_mesh_examples:main',
            'mrc_coarsen = ffeamesh.scripts.mrc_coarsen:main',
            'mrc_crop = ffeamesh.scripts.mrc_crop:main',
            'mrc_header_info = ffeamesh.scripts.mrc_header_info:main',
            'mrc_image_stats = ffeamesh.scripts.mrc_image_stats:main',
            'mrc_to_tets = ffeamesh.scripts.mrc_to_test:main',
            'mrc_threshold = ffeamesh.scripts.mrc_threshold:main',
            'mrc_voxel_size = ffeamesh.scripts.mrc_voxel_size:main',
            'tgv = ffeamesh.app_tgv.tgv:main'
        ]
    }
    )
