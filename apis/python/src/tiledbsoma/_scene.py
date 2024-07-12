# Copyright (c) 2024 TileDB, Inc.
#
# Licensed under the MIT License.

"""Implementation of a SOMA Scene."""

from typing import MutableMapping, Union

from somacore import coordinates, scene

from ._collection import Collection, CollectionBase
from ._coordinates import CompositeTransform, CoordinateSystem
from ._dataframe import DataFrame
from ._dense_nd_array import DenseNDArray
from ._sparse_nd_array import SparseNDArray
from ._soma_object import AnySOMAObject


class Scene(  # type: ignore[misc]  # __eq__ false positive
    CollectionBase[AnySOMAObject],
    scene.Scene[  # type: ignore[type-var]
        DataFrame,
        Collection[
            Union[DenseNDArray, SparseNDArray]
        ],  # not just NDArray since NDArray does not have a common `read`
        AnySOMAObject,
    ],
):
    """TODO: Add documentation for a Scene

    Lifecycle:
        Experimental.
    """

    __slots__ = ()

    _subclass_constrained_soma_types = {
        "img": ("SOMACollection",),
        "obsl": ("SOMACollection",),
        "varl": ("SOMACollection",),
    }

    @property
    def local_coordinate_system(self) -> CoordinateSystem:
        """Coordinate system for this scene."""
        # TODO: Metadata loading should be moved to open/reopen methods.
        coords = self.metadata.get("_soma_local_coordinates")
        if coords is None:
            return CoordinateSystem([])
        coord_sys = CoordinateSystem.from_json(coords)
        return coord_sys

    @property
    def transformations(self) -> MutableMapping[str, coordinates.CoordinateTransform]:
        """Transformations saved for this scene."""
        # TODO: Metadata loading should be moved to open/reopen methods.
        prefix = "_soma_transform_"
        return {
            key[len(prefix) :]: CompositeTransform.from_json(val)
            for key, val in self.metadata.items()
            if key.startswith(prefix)
        }
