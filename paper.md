---
title: 'Tetrahedral meshing for FFEA'
tags:
  - Python
  - video processing
  - crystal growth
  - synchrotron radiation
authors:
  - name: Joanna Leng^[co-first author]
    orcid: 0000-0001-9790-162X
    affiliation: "1"
  - name: Jonathan H. Pickering^[co-first author]
    affiliation: "1"
    orcid: 0000-0001-5283-8065
affiliations:
 - name: School of Computer Science, University of Leeds, Leeds, United Kingdom
   index: 1
date: 31 October 2023
bibliography: paper.bib
---

# Summary

Fluctuating Finite Element Analysis (FFEA) is mesoscale modelling research software package and, in 2020 a combined assessment of software and research interests was performed involving users and developers [?].  This assessment found that meshing was the main bottleneck in the workflow, often taking six months or more. As a result, a Python tetrahedral mesh builder and viewer was developed to suite the mesoscale.

# Statement of Need

FFEA is a continuum mechanics package aimed at modelling the behaviour of proteins and other globular macromolecules [1, 2, 3, 4]. It implements a novel combination of finite element methods and stochastic thermal noise, which raise issues subtly different to tradition finite element methods, some of which relate to the tetrahedral meshes it employs.

A typical system will consist of multiple globular proteins, with attached sugars and potentially water bonded protein surface, these are modelled from cryogenic electron microscope (cryo-em) images in Medical Research Council (MRC) format [5]. Surface refinement techniques used in the structural analyses of manufactured artefacts are not well suited to the irregularities of such systems, and give rise to many small tetrahedra, and very thin tetrahedra, called slither elements. The small and slither elements are problematic as the simulated thermal noise is particularly susceptible to numerical instability, resulting in very short simulation time steps and long calculations, typically 300 000 steps of 5x10-13 seconds. In any case the surfaces of large molecules are not well defined and current surface refinement methods are not ideal.

Originally, the meshes were created from cryo-em images using TetGen [6] then manually edited and refined. Even so the meshes were very unstable, even with very small time steps the simulations crashed after a very short simulation time.  Also, the FFEA research team had problems visualising their meshes, as molecular visualization software does not visualize meshes and 3D visualization software was too complex for the biologists to use and identify the slither elements. A visualization tool was developed to suite these specific needs.

![The viewer’s user interface, showing a meshed model of myosin on the left, with the view controls below. On the right the table of properties of the tetrahedra is can be seen and the global properties of the mesh..\label{fig:example}](images/CGT_drawing_tab.png)

# Design

Since the FFEA project already had several support tools implemented in Python this work was carried out in Python 3.9.  Graphical user interfaces, needed for the mesh viewer, were implemented using PyQt5 and the 3D graphics scene was implemented using pyopengl.  The package was built in a Conda environment, pip install allows useful command line scripts to be provided.

# Functionality

The viewer allows the user to load a tetrahedral mesh from Tetgen or FFEA files. The mesh can be viewed and interactively manipulated in the graphical scene on the left of the application where the user can orient the data to any view. An interactive table to the right displays the properties of the tetrahedra, including shape factor and size, geometrical factors allowing the identification of slither elements. If a tetrahedron is selected in the table it will be highlighted in the viewer.

The most important script converts MRC files to a uniform tetrahedron array using either a five or a six tetrahedron per voxel decomposition. The 5-fold decomposition is not symmetrical so we use the Diamond cubic lattice with alternating left- and right-handed cells. This method of mesh generation not only removes small and slither elements from the final mesh but adds consistency to the meshes.

The command line scripts allow the manipulation of the MRC files that hold the input cryo-em data to:

1. crop a 3D sub-image.
2. print out MRC file header.
3. image intensity statistics and a bin count.
4. copy an MRC file setting all voxels below a threshold to zero.
5. Print out MRC voxel size.

The tetrahedral mesh is constructed from the MRC file, where the user can alter the size of the internal voxels, and threshold out tetrahedra which are not part of the molecule.

The new meshes have been tested by a research student and are under further investigation by the FFEA team. Another research student has explored further optimisations of the surface but with little success in further refining the surface while not reducing the size or shape of the elements significantly.

# For Developers

Developers who wish to contribute to this or to any part of the FFEA software please contact the FFEA team, given in the documentation. The meshing package has been designed as a standalone package, so researchers can use adapt it for their own needs. We recommend researchers who want to extend the package start by copying the repository and setting up the virtual environment.
A test suite is provided. Documentation has been developed using Doxygen to extract documentation from source code comments. The code has been developed using the Pylint static code analysis tool.

# Availability

The software can be obtained from GitHub [GitHub](https://github.com/jonathanHuwP/tet_mesh_tools) under the GNU General Public License (GPL)
Version 3.

# Acknowledgement

Thank you to Molly Garvett, Oliver Harlen, Sarah Harris, Danial Read, Jarvellis Rogers and Yanlai Zhou. This work was supported by the EPSRC grant EP/R025819/1.

# References
