"""SOMA powered by TileDB

SOMA -- stack of matrices, annotated -- is a flexible, extensible, and
open-source API enabling access to data in a variety of formats, and is
motivated by use cases from single-cell biology. The ``tiledbsoma``
Python package is an implementation of SOMA using the
[TileDB Embedded](https://github.com/TileDB-Inc/TileDB) engine.

Provides:
---------
  1. The ability to store, query, and retrieve larger-than-core datasets,
     resident in both cloud (object-store) and local (file) systems.
  2. A data model supporting dataframes, and both sparse and dense
     multi-dimensional arrays.
  3. An extended data model with support for single-cell biology data.

See the [SOMA GitHub repo](https://github.com/single-cell-data/SOMA) for more
information on the SOMA project.

Using the documentation:
------------------------
Coming soon: web-based documentation site.

Documentation is also available via the Python builtin `help` function. We
recommend exploring the package. For example:

>>> import tiledbsoma
>>> help(tiledbsoma.DataFrame)

API maturity tags:
------------------
Classes and functions are annotated with API maturity tags, for example:

    [lifecycle: experimental]

These tags indicate the maturity of each interface, and are patterned after
the RStudio lifecycle stage model. Tags are:
  - experimental: Under active development and may undergo significant and
    breaking changes.
  - maturing: Under active development but the interface and behavior have
    stabilized and are unlikely to change significantly but breaking changes
    are still possible.
  - stable: The interface is considered stable and breaking changes will be
    avoided where possible. Breaking changes that cannot be avoided will be
    accompanied by a major version bump.
  - deprecated: The API is no longer recommended for use and may be removed
    in a future release.

If no tag is present, the assumed state is ``experimental``.

Data types:
-----------
The principle persistent types provided by SOMA are:
  - ``Collection`` - a string-keyed container of SOMA objects.
  - ``DataFrame`` - a multi-column table with a user-defined schema,
    defining the number of columns and their respective column name
    and value type.
  - ``SparseNDArray`` - a sparse multi-dimensional array, storing
    Arrow primitive data types, i.e., int, float, etc.
  - ``DenseNDArray`` -- a dnese multi-dimensional array, storing
    Arrow primitive data types, i.e., int, float, etc.
  - ``Experiment`` - a specialized ``Collection``, representing an
    annotated 2-D matrix of measurements.
  - ``Measurement`` - a specialized ``Collection``, for use within
    the ``Experiment`` class, representing a set of measurements on
    a single set of variables (features, e.g., genes)

SOMA ``Experiment`` and ``Measurement`` are inspired by use cases from
single-cell biology.

SOMA uses the [Arrow](https://arrow.apache.org/docs/python/index.html) type
system and memory model for its in-memory type system and schema. For
example, the schema of a ``tiledbsoma.DataFrame`` is expressed as an
[Arrow Schema](https://arrow.apache.org/docs/python/data.html#schemas).
"""

from somacore import AxisColumnNames, AxisQuery, ExperimentAxisQuery
from somacore.options import ResultOrder

from ._collection import Collection
from ._constants import SOMA_JOINID
from ._dataframe import DataFrame
from ._dense_nd_array import DenseNDArray
from ._exception import DoesNotExistError, SOMAError
from ._experiment import Experiment
from ._factory import open
from ._general_utilities import (
    get_implementation,
    get_implementation_version,
    get_SOMA_version,
    get_storage_engine,
    show_package_versions,
)
from ._measurement import Measurement
from ._sparse_nd_array import SparseNDArray
from .libtiledbsoma import (
    tiledbsoma_stats_disable,
    tiledbsoma_stats_dump,
    tiledbsoma_stats_enable,
    tiledbsoma_stats_reset,
)
from .options import SOMATileDBContext

__version__ = get_implementation_version()

__all__ = [
    "AxisColumnNames",
    "AxisQuery",
    "Collection",
    "DataFrame",
    "DenseNDArray",
    "DoesNotExistError",
    "Experiment",
    "ExperimentAxisQuery",
    "get_implementation_version",
    "get_implementation",
    "get_SOMA_version",
    "get_storage_engine",
    "Measurement",
    "open",
    "show_package_versions",
    "SOMAError",
    "ResultOrder",
    "SOMATileDBContext",
    "SOMA_JOINID",
    "SparseNDArray",
    "tiledbsoma_stats_disable",
    "tiledbsoma_stats_dump",
    "tiledbsoma_stats_enable",
    "tiledbsoma_stats_reset",
]
