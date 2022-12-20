#!/usr/bin/env python
"""
setup.py

File to compile c and prepare python for a pip install.

---------------------------------------


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
from setuptools import setup, find_packages

DISTUTILS_DEBUG=1

setup(
    name='ffeamesh',
    version='1.0.0',
    author="The FFEA Team",
    description ='Meshing tools for the Fluctuating Finite Element Analysis tool (FFEA).',
    url='http://ffea.bitbucket.com',
    license='GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007',
    packages=['ffeamesh'],
    #packages=find_packages(),
    install_requires=[
        'argparse',
        'numpy',
        'vtk',
        'mrcfile'
    ],
    scripts=['ffeamesh/scripts/five_tets.py',
             'ffeamesh/scripts/six_tets.py',
             'ffeamesh/scripts/six_tets_demo.py',
             'ffeamesh/scripts/zoom.py',
             'ffeamesh/scripts/fft_smooth.py',
             'ffeamesh/scripts/make_test_mrcfile.py',
             'ffeamesh/scripts/mrc_threshold.py',
             'ffeamesh/scripts/mrc_crop.py',
             'ffeamesh/scripts/mrc_header_info.py',
             'ffeamesh/scripts/mrc_image_stats.py',
             'tests/tetvolume.py'])
