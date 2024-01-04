#!/usr/bin/env python
"""
pip install the scripts

You should have received a copy of the GNU General Public License
along with FFEA.  If not, see <http://www.gnu.org/licenses/>.

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
            'mrc_to_tets = ffeamesh.scripts.mrc_to_tets:main',
            'mrc_threshold = ffeamesh.scripts.mrc_threshold:main',
            'mrc_voxel_size = ffeamesh.scripts.mrc_voxel_size:main',
            'tgv = ffeamesh.app_tgv.tgv:main'
        ]
    }
    )
