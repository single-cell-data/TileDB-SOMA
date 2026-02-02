# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.
"""Implementations of the composed SOMA data types."""

from typing import Final, Generic, TypeVar

from . import _mixin, base, collection, data

_DF = TypeVar("_DF", bound=data.DataFrame)
"""A particular implementation of DataFrame."""
_NDColl = TypeVar("_NDColl", bound=collection.Collection[data.NDArray])
"""A particular implementation of a collection of NDArrays."""
_DenseNDColl = TypeVar("_DenseNDColl", bound=collection.Collection[data.DenseNDArray])
"""A particular implementation of a collection of DenseNDArrays."""
_SparseNDColl = TypeVar("_SparseNDColl", bound=collection.Collection[data.SparseNDArray])
"""A particular implementation of a collection of SparseNDArrays."""
_RootSO = TypeVar("_RootSO", bound=base.SOMAObject)
"""The root SomaObject type of the implementation."""


class Measurement(
    collection.BaseCollection[_RootSO],
    Generic[_DF, _NDColl, _DenseNDColl, _SparseNDColl, _RootSO],
):
    """A set of observations defined by a dataframe, with measurements.

    This is a common set of annotated variables (defined by the ``var``
    dataframe) for which values (e.g., measurements or calculations) are stored
    in sparse and dense ND arrays.

    The observables are inherited from the parent ``Experiment``'s
    ``obs`` dataframe. The ``soma_joinid`` of these observables (``obsid``),
    along with those of the measurement's ``var`` dataframe (``varid``),
    are the indices for all the other matrices stored in the measurement.

    Lifecycle: maturing
    """

    __slots__ = ()
    soma_type: Final = "SOMAMeasurement"  # type: ignore[misc]

    var = _mixin.item[_DF]()
    """Primary annotations on the variable axis for vars on this measurement.

    This annotates _columns_ of the ``X`` arrays. The contents of the
    ``soma_joinid`` pseudo-column define the variable index domain (``varid``)
    All variables for this measurement _must_ be defined in this dataframe.
    """

    X = _mixin.item[_NDColl]()
    """A collection of matrices containing feature values.

    Each matrix is indexed by ``[obsid, varid]``. Sparse and dense 2D arrays may
    both be used in any combination in ``X``.
    """

    obsm = _mixin.item[_NDColl]()
    """Matrices containing annotations of each ``obs`` row.

    This has the same shape as ``obs`` and is indexed with ``obsid``.
    """

    obsp = _mixin.item[_SparseNDColl]()
    """Matrices containing pairwise annotations of each ``obs`` row.

    This is indexed by ``[obsid_1, obsid_2]``.
    """

    varm = _mixin.item[_NDColl]()
    """Matrices containing annotations of each ``var`` row.

    This has the same shape as ``var`` and is indexed with ``varid``.
    """

    varp = _mixin.item[_SparseNDColl]()
    """Matrices containing pairwise annotations of each ``var`` row.

    This is indexed by ``[varid_1, varid_2]``.
    """

    var_spatial_presence = _mixin.item[_DF]()
    """A dataframe that stores the presence of var in the spatial scenes.

    This provides a join table for the var ``soma_joinid`` and the scene names used in
    the ``spatial`` collection. This dataframe must contain index columns ``soma_joinid``
    and ``scene_id``. The ``scene_id`` column  must have type ``string``. The
    dataframe must contain a ``boolean`` column ``data``. The values of ``data`` are
    ``True`` if the var with varid ``soma_joinid`` is contained in scene with name
    ``scene_id`` and ``False`` otherwise.
    """
