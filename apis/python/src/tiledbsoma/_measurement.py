# Copyright (c) 2021-2023 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2023 TileDB, Inc.
#
# Licensed under the MIT License.

"""Implementation of a SOMA Measurement.
"""

from typing import Union

from somacore import measurement

from ._collection import Collection, CollectionBase
from ._dataframe import DataFrame
from ._dense_nd_array import DenseNDArray
from ._sparse_nd_array import SparseNDArray
from ._tiledb_object import AnyTileDBObject


class Measurement(  # type: ignore[misc]  # __eq__ false positive
    CollectionBase[AnyTileDBObject],
    measurement.Measurement[  # type: ignore[type-var]
        DataFrame,
        Collection[
            Union[SparseNDArray, DenseNDArray]
        ],  # not just `NDArray` since that has no common `read`
        Collection[DenseNDArray],
        Collection[SparseNDArray],
        AnyTileDBObject,
    ],
):
    """A set of observations defined by a dataframe, with measurements.

    This is a common set of annotated variables (defined by the ``var``
    dataframe) for which values (e.g., measurements or calculations) are stored
    in sparse and dense ND arrays.

    The observables are inherited from the parent ``Experiment``'s
    ``obs`` dataframe. The ``soma_joinid`` of these observables (``obsid``),
    along with those of the measurement's ``var`` dataframe (``varid``),
    are the indices for all the other matrices stored in the measurement.

    In most cases, users interact with a measurement via querying the experiment
    which contains it, rather than directly accessing its fields.

    Attributes:
        var (DataFrame):
            Primary annotations on the variable axis, for variables in this measurement
            (i.e., annotates columns of ``X``). The contents of the ``soma_joinid`` column
            define the variable index domain, AKA var_id. All variables for this measurement
            must be defined in this dataframe.
        X (Collection[SparseNDArray]):
            A collection of sparse matrices, each containing measured feature values.
            Each matrix is indexed by ``[obsid, varid]``.
        obsm (Collection[DenseNDArray]):
            A collection of dense matrices containing annotations of each ``obs`` row.
            Has the same shape as ``obs``, and is indexed with ``obsid``.
        obsp (Collection[SparseNDArray]):
            A collection of sparse matrices containing pairwise annotations of each ``obs`` row.
            Indexed with ``[obsid_1, obsid_2]``.
        varm (Collection[DenseNDArray]):
            A collection of dense matrices containing annotations of each ``var`` row.
            Has the same shape as ``var``, and is indexed with ``varid``.
        varp (Collection[SparseNDArray]):
            A collection of sparse matrices containing pairwise annotations of each ``var`` row.
            Indexed with ``[varid_1, varid_2]``

    Lifecycle:
        Experimental.
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
