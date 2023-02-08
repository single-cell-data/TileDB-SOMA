from somacore import measurement

from .collection import Collection, CollectionBase
from .common_nd_array import NDArray
from .dataframe import DataFrame
from .dense_nd_array import DenseNDArray
from .sparse_nd_array import SparseNDArray
from .tiledb_object import AnyTileDBObject


class Measurement(
    CollectionBase[AnyTileDBObject],
    measurement.Measurement[  # type: ignore[type-var]
        DataFrame,
        Collection[NDArray],
        Collection[DenseNDArray],
        Collection[SparseNDArray],
        AnyTileDBObject,
    ],
):
    """
    A ``Measurement`` is a sub-element of a ``Experiment``, and is otherwise
    a specialized ``Collection`` with pre-defined fields:

    ``var``: ``DataFrame``

    Primary annotations on the variable axis, for variables in this measurement
    (i.e., annotates columns of ``X``). The contents of the ``soma_joinid`` column
    define the variable index domain, AKA var_id. All variables for this measurement
    must be defined in this dataframe.

    ``X``: ``Collection`` of ``SparseNDArray``

    A collection of sparse matrices, each containing measured feature values.
    Each matrix is indexed by ``[obsid, varid]``.

    ``obsm``: ``Collection`` of ``DenseNDArray``

    A collection of dense matrices containing annotations of each ``obs`` row.
    Has the same shape as ``obs``, and is indexed with ``obsid``.

    ``obsp``: ``Collection`` of ``SparseNDArray``

    A collection of sparse matrices containing pairwise annotations of each ``obs`` row.
    Indexed with ``[obsid_1, obsid_2]``.

    ``varm``: ``Collection`` of ``DenseNDArray``

    A collection of dense matrices containing annotations of each ``var`` row.
    Has the same shape as ``var``, and is indexed with ``varid``.

    ``varp``: ``Collection`` of ``SparseNDArray``

    A collection of sparse matrices containing pairwise annotations of each ``var`` row.
    Indexed with ``[varid_1, varid_2]``

    [lifecycle: experimental]
    """

    __slots__ = ()

    _subclass_constrained_soma_types = {
        "var": ("SOMADataFrame",),
        "X": ("SOMACollection",),
        "obsm": ("SOMACollection",),
        "obsp": ("SOMACollection",),
        "varm": ("SOMACollection",),
        "varp": ("SOMACollection",),
    }
