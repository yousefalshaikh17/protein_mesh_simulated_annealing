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

 @author: jonathan pickering, joanna leng 26 Jan 23
"""
import tetmeshtools.meshtools.tetgenstructs as ts
import tetmeshtools.meshtools.tetgenread as tr

def test_uniformity():
    """
    test uniformity
    """
    values = [2, 3, 4]
    print(f"Values {values} => {tr.uniformity(values)} should be 0.07407")

def test_edges_to_area_ratio():
    """
    test edges_to_area_ratio
    """
    nodes = []
    nodes.append(ts.NodePoint(7,   0.0, 0.0, 0.0))
    nodes.append(ts.NodePoint(12,  1.0, 0.0, 0.0))
    nodes.append(ts.NodePoint(123, 0.0, 1.0, 0.0))
    print(f"Edges to area ratio {tr.edges_to_area_ratio(nodes)} should be 23.3137")
