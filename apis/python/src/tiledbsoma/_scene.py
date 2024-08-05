# Copyright (c) 2024 TileDB, Inc.
#
# Licensed under the MIT License.

"""Implementation of a SOMA Scene."""


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

    __slots__ = ()
    _wrapper_type = _tdb_handles.SceneWrapper

    _subclass_constrained_soma_types = {
        "img": ("SOMACollection",),
        "obsl": ("SOMACollection",),
        "varl": ("SOMACollection",),
    }

    @property
    def coordinate_space(self) -> CoordinateSpace:
        """Coordinate system for this scene."""
        # TODO: Metadata loading should be moved to open/reopen methods.
        coords = self.metadata.get("_soma_local_coordinates")
        if coords is None:
            return CoordinateSpace([])
        return CoordinateSpace.from_json(coords)
