# ----------------------------------------------------------------
# The way this works is:
#
# o We put a `[tool.setuptools_scm]` stanza in pyproject.toml.
#
# o When the user does `pip install .` the src/tiledbsoma/_version.py
#   file gets created. We have a .gitignore line for that since it
#   isn't supposed to go into version control, but rather, is supposed
#   to be created at package-setup time.
#
# o The src/tiledbsoma/_version.py is created by setuptools by reading
#   Git history. If the currently checked-out code is at a version stamp
#   like `1.2.3` we'll get `__version__` being `1.2.3`; else we'll get
#   something like `1.2.4.dev2`.
#
# o When the user does `import tiledbsoma` this code here is run,
#   which loads the `src/tiledbsoma/_version.py`.
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
        __version__ = get_distribution("tiledbsoma").version
    except DistributionNotFound:
        __version__ = "unknown"
# ----------------------------------------------------------------

from .collection import Collection
from .dataframe import DataFrame
from .dense_nd_array import DenseNdArray
from .exception import DoesNotExistError, Error
from .experiment import Experiment
from .general_utilities import (
    get_implementation,
    get_implementation_version,
    get_SOMA_version,
    get_storage_engine,
)
from .indexed_dataframe import IndexedDataFrame
from .measurement import Measurement
from .metadata_mapping import MetadataMapping
from .query_condition import QueryCondition  # type: ignore
from .sparse_nd_array import SparseNdArray
from .tiledb_array import TileDBArray
from .tiledb_object import TileDBObject
from .tiledb_platform_config import TileDBPlatformConfig

__all__ = [
    "get_implementation",
    "get_implementation_version",
    "get_SOMA_version",
    "get_storage_engine",
    "TileDBPlatformConfig",
    "TileDBObject",
    "TileDBArray",
    "Collection",
    "DataFrame",
    "DenseNdArray",
    "DoesNotExistError",
    "Experiment",
    "Error",
    "IndexedDataFrame",
    "Measurement",
    "MetadataMapping",
    "SparseNdArray",
    "QueryCondition",
]
