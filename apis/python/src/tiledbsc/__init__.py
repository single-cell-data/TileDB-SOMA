try:
    from ._version import version as __version__
except ImportError:
    from pkg_resources import DistributionNotFound, get_distribution
    try:
        __version__ = get_distribution("tiledbsc").version
    except DistributionNotFound:
        __version__ = "unknown"

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
