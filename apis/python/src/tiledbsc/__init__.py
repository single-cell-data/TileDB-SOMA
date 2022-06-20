# This reads ../../pyproject.toml
import importlib.metadata

__version__ = importlib.metadata.version("tiledbsc")

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
