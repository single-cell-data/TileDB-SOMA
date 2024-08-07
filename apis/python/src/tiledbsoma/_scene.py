# Copyright (c) 2024 TileDB, Inc.
#
# Licensed under the MIT License.

"""Implementation of a SOMA Scene."""

from typing import Any, Optional

from somacore import scene

from . import _tdb_handles
from ._collection import Collection, CollectionBase
from ._coordinates import (
    Axis,
    CoordinateSpace,
    CoordinateTransform,
    IdentityCoordinateTransform,
    transform_to_json,
)
from ._exception import SOMAError
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
            self._coord_space: Optional[CoordinateSpace] = None
        else:
            self._coord_space = CoordinateSpace.from_json(coord_space)

    @property
    def coordinate_space(self) -> Optional[CoordinateSpace]:
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
        self._coord_space = None

    def register_image2d(
        self,
        image_name: str,
        transform: CoordinateTransform,
        *,
        subcollection_name: str = "img",
        coordinate_space: Optional[CoordinateSpace] = None,
    ) -> Image2DCollection:
        """TODO Add docstring"""
        # Check the transform matches this
        if self.coordinate_space is None:
            raise SOMAError(
                "The scene coordinate space must be set before registering an image."
            )
        if transform.output_axes != self.coordinate_space.axis_names:
            raise ValueError(
                f"The name of the transform output axes, {transform.output_axes}, do "
                f"not match the name of the axes in the scene coordinate space, "
                f"{self.coordinate_space.axis_names}."
            )

        # Create the coordinate space if it does not exist. Otherwise, check it is
        # compatible with the provide transform.
        if coordinate_space is None:
            if isinstance(transform, IdentityCoordinateTransform):
                coordinate_space = self.coordinate_space
            else:
                coordinate_space = CoordinateSpace(
                    tuple(Axis(name=axis_name) for axis_name in transform.input_axes)
                )
        else:
            if transform.input_axes != coordinate_space.axis_names:
                raise ValueError(
                    f"The name of the transform input axes, {transform.input_axes}, do "
                    f"not match the name of the axes in the provided coordinate space, "
                    f"{coordinate_space.axis_names}."
                )

        # Check asset exists in the specified location.
        try:
            subcollection: Collection = self[subcollection_name]  # type: ignore
        except KeyError as ke:
            raise KeyError(
                f"No collection '{subcollection_name}' in this scene."
            ) from ke
        try:
            image: Image2DCollection = subcollection[image_name]
        except KeyError as ke:
            raise KeyError(
                f"No Image2DCollection named '{image_name}' in '{subcollection_name}'."
            ) from ke
        if not isinstance(image, Image2DCollection):
            raise TypeError(
                f"'{image_name}' in '{subcollection_name}' is not an Image2DCollection."
            )

        image.coordinate_space = coordinate_space
        subcollection.metadata[f"soma_scene_registry_{image_name}"] = transform_to_json(
            transform
        )
        return image
