# ----------------------------------------------------------------
import setuptools_scm

# __version__ = setuptools_scm.get_version(root='../../../..')
# Nope:
# >>> import tiledbsc
# >>> tiledbsc.__version__
# LookupError: setuptools-scm was unable to detect version for /Users/johnkerl/git.

# __version__ = setuptools_scm.get_version(root='../../..')
# >>> import tiledbsc
# >>> tiledbsc.__version__
# Nope:
# LookupError: setuptools-scm was unable to detect version for /Users/johnkerl/git/single-cell-data.

__version__ = setuptools_scm.get_version(root='../..')
# >>> import tiledbsc
# >>> tiledbsc.__version__
# '0.0.2.dev0+g1e1e691.d20220620'
# BUT this only works when I'm cd'ed into /Users/johnkerl/git/single-cell-data/TileDB-SingleCell/apis/python

# ----------------------------------------------------------------
from .soma_collection import SOMACollection
from .soma import SOMA
from .soma_options import SOMAOptions

from .tiledb_object import TileDBObject
from .tiledb_array import TileDBArray
from .tiledb_group import TileDBGroup
from .assay_matrix import AssayMatrix
from .annotation_matrix import AnnotationMatrix

from .annotation_matrix_group import AnnotationMatrixGroup
from .annotation_pairwise_matrix_group import AnnotationPairwiseMatrixGroup
from .assay_matrix_group import AssayMatrixGroup
from .raw_group import RawGroup
from .uns_group import UnsGroup
from .uns_array import UnsArray

from .util_ann import describe_ann_file
from .util_tiledb import show_single_cell_group
