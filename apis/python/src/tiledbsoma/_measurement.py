# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""Implementation of a SOMA Measurement."""

from typing import Final, Union

from . import _mixin
from . import pytiledbsoma as clib
from ._collection import Collection, CollectionBase
from ._dataframe import DataFrame
from ._dense_nd_array import DenseNDArray
from ._soma_object import SOMAObject
from ._sparse_nd_array import SparseNDArray


class Measurement(CollectionBase[SOMAObject]):
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
        Maturing.
    """

    __slots__ = ()
    _handle_type = clib.SOMAMeasurement
    soma_type: Final = "SOMAMeasurement"  # type: ignore[misc]

    var = _mixin.item[DataFrame]()
    """Primary annotations on the variable axis for vars on this measurement.

    This annotates _columns_ of the ``X`` arrays. The contents of the
    ``soma_joinid`` pseudo-column define the variable index domain (``varid``)
    All variables for this measurement _must_ be defined in this dataframe.
    """

    X = _mixin.item[Collection[Union[SparseNDArray, DenseNDArray]]]()
    """A collection of matrices containing feature values.

    Each matrix is indexed by ``[obsid, varid]``. Sparse and dense 2D arrays may
    both be used in any combination in ``X``.
    """

    obsm = _mixin.item[Collection[Union[SparseNDArray, DenseNDArray]]]()
    """Matrices containing annotations of each ``obs`` row.

    This has the same shape as ``obs`` and is indexed with ``obsid``.
    """

    obsp = _mixin.item[Collection[SparseNDArray]]()
    """Matrices containing pairwise annotations of each ``obs`` row.

    This is indexed by ``[obsid_1, obsid_2]``.
    """

    varm = _mixin.item[Collection[Union[SparseNDArray, DenseNDArray]]]()
    """Matrices containing annotations of each ``var`` row.

    This has the same shape as ``var`` and is indexed with ``varid``.
    """

    varp = _mixin.item[Collection[SparseNDArray]]()
    """Matrices containing pairwise annotations of each ``var`` row.

    This is indexed by ``[varid_1, varid_2]``.
    """

    var_spatial_presence = _mixin.item[DataFrame]()
    """A dataframe that stores the presence of var in the spatial scenes.

    This provides a join table for the var ``soma_joinid`` and the scene names used in
    the ``spatial`` collection. This dataframe must contain index columns ``soma_joinid``
    and ``scene_id``. The ``scene_id`` column  must have type ``string``. The
    dataframe must contain a ``boolean`` column ``data``. The values of ``data`` are
    ``True`` if the var with varid ``soma_joinid`` is contained in scene with name
    ``scene_id`` and ``False`` otherwise.
    """
