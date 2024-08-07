# Copyright (c) 2024 TileDB, Inc.
#
# Licensed under the MIT License.

"""Implementation of a SOMA Scene."""

from typing import Any

from somacore import scene

from . import _tdb_handles
from ._collection import Collection, CollectionBase
from ._coordinates import CoordinateSpace
from ._images import Image2DCollection
from ._soma_object import AnySOMAObject
from ._spatial_dataframe import SpatialDataFrame


class Scene(  # type: ignore[misc]  # __eq__ false positive
    CollectionBase[AnySOMAObject],
    scene.Scene[  # type: ignore[type-var]
        Collection[SpatialDataFrame], Image2DCollection, AnySOMAObject
    ],
):
    """TODO: Add documentation for a Scene

    Lifecycle:
        Experimental.
    """

    __slots__ = "_coord_space"
    _wrapper_type = _tdb_handles.SceneWrapper

    _subclass_constrained_soma_types = {
        "img": ("SOMACollection",),
        "obsl": ("SOMACollection",),
        "varl": ("SOMACollection",),
    }

    def __init__(
        self,
        handle: _tdb_handles.SOMAGroupWrapper[Any],
        **kwargs: Any,
    ):
        super().__init__(handle, **kwargs)
        coord_space = self.metadata.get("soma_coordinate_space")
        if coord_space is None:
            self._coord_space = CoordinateSpace([])
        else:
            self._coord_space = CoordinateSpace.from_json(coord_space)

    @property
    def coordinate_space(self) -> CoordinateSpace:
        """Coordinate system for this scene."""
        return self._coord_space

    @coordinate_space.setter
    def coordinate_space(self, value: CoordinateSpace) -> None:
        if not isinstance(value, CoordinateSpace):
            raise TypeError(f"Invalid type {type(value).__name__}.")
        self.metadata["soma_coordinate_space"] = value.to_json()
        self._coord_space = value

    @coordinate_space.deleter
    def coordinate_space(self) -> None:
        del self.metadata["soma_coordinate_space"]
        self._coord_space = CoordinateSpace([])
