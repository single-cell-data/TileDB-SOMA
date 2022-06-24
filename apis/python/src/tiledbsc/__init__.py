# The way this works is:
#
# o We put a `[tool.setuptools_scm]` stanza in pyproject.toml.
#
# o When the user does `pip install .` the src/tiledbsc/_version.py
#   file gets created. We have a .gitignore line for that since it
#   isn't supposed to go into version control, but rather, is supposed
#   to be created at package-setup time.
#
# o The src/tiledbsc/_version.py is created by setuptools by reading
#   Git history. If the currently checked-out code is at a version stamp
#   like `1.2.3` we'll get `__version__` being `1.2.3`; else we'll get
#   something like `1.2.4.dev2`.
#
# o When the user does `import tiledbsc` this code here is run,
#   which loads the `src/tiledbsc/_version.py`.
#
# In summary, who needs to do what:
#
# o Package developers: simply create a Git release at
# https://github.com/single-cell-data/TileDB-SingleCell/releases
#
# o Package users: just `pip install .`.
try:
    from ._version import version as __version__
except ImportError:
    from pkg_resources import DistributionNotFound, get_distribution

    try:
        __version__ = get_distribution("tiledbsc").version
    except DistributionNotFound:
        __version__ = "unknown"

import logging

from .annotation_matrix import AnnotationMatrix
from .annotation_matrix_group import AnnotationMatrixGroup
from .annotation_pairwise_matrix_group import AnnotationPairwiseMatrixGroup
from .assay_matrix import AssayMatrix
from .assay_matrix_group import AssayMatrixGroup
from .raw_group import RawGroup
from .soma import SOMA
from .soma_collection import SOMACollection
from .soma_options import SOMAOptions
from .soma_slice import SOMASlice
from .tiledb_array import TileDBArray
from .tiledb_group import TileDBGroup
from .tiledb_object import TileDBObject
from .uns_array import UnsArray
from .uns_group import UnsGroup
from .util_ann import describe_ann_file
from .util_tiledb import show_soma_schemas
