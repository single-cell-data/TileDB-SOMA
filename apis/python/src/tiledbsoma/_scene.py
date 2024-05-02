# Copyright (c) 2024 TileDB, Inc.
#
# Licensed under the MIT License.

"""Implementation of a SOMA Scene."""


from typing import Union

from somacore import scene

from ._collection import Collection, CollectionBase
from ._dataframe import DataFrame
from ._dense_nd_array import DenseNDArray
from ._sparse_nd_array import SparseNDArray
from ._tiledb_object import AnyTileDBObject


class Scene(  # type: ignore[misc]  # __eq__ false positive
    CollectionBase[AnyTileDBObject],
    scene.Scene[  # type: ignore[type-var]
        Collection[
            Union[DataFrame, DenseNDArray, SparseNDArray]
        ],  # not just DataFrame and NDArray since NDArray does not have a common `read`
        AnyTileDBObject,
    ],
):
    """TODO: Add documentation for a Scene

    Lifecycle:
        Experimental.
    """

    __slots__ = ()

    _subclass_constrained_soma_types = {
        "exp": ("SOMACollection",),
        "ms": ("SOMACollection",),
    }
