#!/usr/bin/env python
"""
pip install the scripts

You should have received a copy of the GNU General Public License.
If not, see <http://www.gnu.org/licenses/>.

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

This work was funded by Joanna Leng's EPSRC funded RSE Fellowship (EP/R025819/1)

@copyright 2020
@author: j.h.pickering@leeds.ac.uk and j.leng@leeds.ac.uk
"""
from setuptools import setup

DISTUTILS_DEBUG=1

setup(
    name='tetmeshtools',
    version='1.0.0',
    author="J. Leng, J. Pickering",
    description ='Meshing tools developed for the Fluctuating Finite Element Analysis tool (FFEA).',
    url='https://github.com/jonathanHuwP/tet_mesh_tools',
    license='GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007',
    packages=['tetmeshtools'],
    install_requires=[
        'argparse',
        'numpy',
        'vtk',
        'mrcfile'
    ],
    entry_points={
        'console_scripts': [
            'make_test_mrc_file = tetmeshtools.scripts.make_test_mrcfile:main',
            'make_tet_mesh_examples = tetmeshtools.scripts.make_tet_mesh_examples:main',
            'mrc_coarsen = tetmeshtools.scripts.mrc_coarsen:main',
            'mrc_crop = tetmeshtools.scripts.mrc_crop:main',
            'mrc_header_info = tetmeshtools.scripts.mrc_header_info:main',
            'mrc_image_stats = tetmeshtools.scripts.mrc_image_stats:main',
            'mrc_to_tets = tetmeshtools.scripts.mrc_to_tets:main',
            'mrc_threshold = tetmeshtools.scripts.mrc_threshold:main',
            'mrc_voxel_size = tetmeshtools.scripts.mrc_voxel_size:main',
            'pdb_to_mrc = tetmeshtools.scripts.pdb_to_mrc:main',
            'tgv = tetmeshtools.app_tgv.tgv:main'
        ]
    }
    )
