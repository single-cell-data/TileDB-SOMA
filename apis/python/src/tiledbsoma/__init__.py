"""SOMA powered by TileDB

SOMA --- stack of matrices, annotated --- is a flexible, extensible, and
open-source API enabling access to data in a variety of formats, and is
motivated by use cases from single-cell biology. The ``tiledbsoma``
Python package is an implementation of SOMA using the
`TileDB Embedded <https://github.com/TileDB-Inc/TileDB>`_ engine.

Provides
----------
  1. The ability to store, query, and retrieve larger-than-core datasets,
     resident in both cloud (object-store) and local (file) systems.
  2. A data model supporting dataframes, and both sparse and dense
     multi-dimensional arrays.
  3. An extended data model with support for single-cell biology data.

See the `SOMA GitHub repo <https://github.com/single-cell-data/SOMA>`_ for more
information on the SOMA project.

Using the documentation
-------------------------

Documentation is also available via the Python builtin ``help`` function. We
recommend exploring the package. For example:

>>> import tiledbsoma
>>> help(tiledbsoma.DataFrame)

API maturity tags
------------------

Classes and functions are annotated with API maturity tags, for example:

    ``Lifecycle: experimental``

These tags indicate the maturity of each interface, and are patterned after
the RStudio lifecycle stage model. Tags are:

- ``experimental``: Under active development and may undergo significant and
  breaking changes.
- ``maturing``: Under active development but the interface and behavior have
  stabilized and are unlikely to change significantly but breaking changes
  are still possible.
- ``stable``: The interface is considered stable and breaking changes will be
  avoided where possible. Breaking changes that cannot be avoided will be
  accompanied by a major version bump.
- ``deprecated``: The API is no longer recommended for use and may be removed
  in a future release.

If no tag is present, the state is ``experimental``.

Data types
------------

The principal persistent types provided by SOMA are:

- :class:`Collection` -- a string-keyed container of SOMA objects.
- :class:`DataFrame` -- a multi-column table with a user-defined schema,
  defining the number of columns and their respective column name
  and value type.
- :class:`SparseNDArray` -- a sparse multi-dimensional array, storing
  Arrow primitive data types, i.e., int, float, etc.
- :class:`DenseNDArray` -- a dense multi-dimensional array, storing
  Arrow primitive data types, i.e., int, float, etc.
- :class:`Experiment` -- a specialized :class:`Collection`, representing an
  annotated 2-D matrix of measurements.
- :class:`Measurement` -- a specialized :class:`Collection`, for use within
  the :class:`Experiment` class, representing a set of measurements on
  a single set of variables (features, e.g., genes)

SOMA :class:`Experiment` and :class:`Measurement` are inspired by use cases from
single-cell biology.

SOMA uses the `Arrow <https://arrow.apache.org/docs/python/index.html>`_ type
system and memory model for its in-memory type system and schema. For
example, the schema of a :class:`DataFrame` is expressed as an
`Arrow Schema <https://arrow.apache.org/docs/python/generated/pyarrow.Schema.html>`_.

Error handling
---------------
Most errors will be signaled with a raised Exception. Of note:

- :class:`NotImplementedError` will be raised when the requested function or method
  is unsupported.
- :class:`SOMAError` is a base class for all SOMA-specific errors.
- ``TileDBError`` will be raised for many TileDB-specific errors.

Most errors will raise an appropriate Python error, e.g., ::class:`TypeError` or
:class:`ValueError`.
"""

# ^^ the rest is autogen whether viewed from Python on-line help, Sphinx/readthedocs, etc.  It's
# crucial that we include a separator (e.g. "Classes and functions") to make an entry in the
# readthedocs table of contents.

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
from .options import SOMATileDBContext, TileDBCreateOptions
from .pytiledbsoma import (
    tiledbsoma_stats_disable,
    tiledbsoma_stats_dump,
    tiledbsoma_stats_enable,
    tiledbsoma_stats_reset,
)

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
    "ResultOrder",
    "show_package_versions",
    "SOMA_JOINID",
    "SOMAError",
    "SOMATileDBContext",
    "SparseNDArray",
    "TileDBCreateOptions",
    "tiledbsoma_stats_disable",
    "tiledbsoma_stats_dump",
    "tiledbsoma_stats_enable",
    "tiledbsoma_stats_reset",
]
