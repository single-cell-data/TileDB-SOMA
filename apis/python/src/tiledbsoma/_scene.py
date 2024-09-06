# Copyright (c) 2024 TileDB, Inc.
#
# Licensed under the MIT License.

"""Implementation of a SOMA Scene."""

from typing import Any, Optional, Union

from somacore import (
    Axis,
    CoordinateSpace,
    CoordinateTransform,
    IdentityTransform,
    scene,
)

from . import _tdb_handles
from ._collection import Collection, CollectionBase
from ._constants import SOMA_COORDINATE_SPACE_METADATA_KEY
from ._coordinates import (
    coordinate_space_from_json,
    coordinate_space_to_json,
    transform_from_json,
    transform_to_json,
)
from ._exception import SOMAError
from ._geometry_dataframe import GeometryDataFrame
from ._multiscale_image import MultiscaleImage
from ._point_cloud import PointCloud
from ._soma_object import AnySOMAObject


class Scene(  # type: ignore[misc]  # __eq__ false positive
    CollectionBase[AnySOMAObject],
    scene.Scene[  # type: ignore[type-var]
        MultiscaleImage,
        Collection[Union[PointCloud, GeometryDataFrame]],
        Collection[Collection[Union[PointCloud, GeometryDataFrame]]],
        AnySOMAObject,
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
        coord_space = self.metadata.get(SOMA_COORDINATE_SPACE_METADATA_KEY)
        if coord_space is None:
            self._coord_space: Optional[CoordinateSpace] = None
        else:
            self._coord_space = coordinate_space_from_json(coord_space)

    @property
    def coordinate_space(self) -> Optional[CoordinateSpace]:
        """Coordinate system for this scene."""
        return self._coord_space

    @coordinate_space.setter
    def coordinate_space(self, value: CoordinateSpace) -> None:
        if not isinstance(value, CoordinateSpace):
            raise TypeError(f"Invalid type {type(value).__name__}.")
        self.metadata[SOMA_COORDINATE_SPACE_METADATA_KEY] = coordinate_space_to_json(
            value
        )
        self._coord_space = value

    def register_point_cloud(
        self,
        point_cloud_name: str,
        transform: CoordinateTransform,
        *,
        subcollection_name: str = "obsl",
        coordinate_space: Optional[CoordinateSpace] = None,
    ) -> PointCloud:
        if self.coordinate_space is None:
            raise SOMAError(
                "The scene coordinate space must be set before registering a point "
                "cloud."
            )
        # Create the coordinate space if it does not exist. Otherwise, check it is
        # compatible with the provide transform.
        if coordinate_space is None:
            if isinstance(transform, IdentityTransform):
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
            point_cloud: PointCloud = subcollection[point_cloud_name]
        except KeyError as ke:
            raise KeyError(
                f"No PointCloud named '{point_cloud_name}' in '{subcollection_name}'."
            ) from ke
        if not isinstance(point_cloud, PointCloud):
            raise TypeError(
                f"'{point_cloud_name}' in '{subcollection_name}' is not an PointCloud."
            )

        point_cloud.coordinate_space = coordinate_space
        subcollection.metadata[f"soma_scene_registry_{point_cloud_name}"] = (
            transform_to_json(transform)
        )
        return point_cloud

    def register_multiscale_image(
        self,
        image_name: str,
        transform: CoordinateTransform,
        *,
        subcollection_name: str = "img",
        coordinate_space: Optional[CoordinateSpace] = None,
    ) -> MultiscaleImage:
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
            if isinstance(transform, IdentityTransform):
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
            image: MultiscaleImage = subcollection[image_name]
        except KeyError as ke:
            raise KeyError(
                f"No multiscale image named '{image_name}' in '{subcollection_name}'."
            ) from ke
        if not isinstance(image, MultiscaleImage):
            raise TypeError(
                f"'{image_name}' in '{subcollection_name}' is not an MultiscaleImage."
            )

        image.coordinate_space = coordinate_space
        subcollection.metadata[f"soma_scene_registry_{image_name}"] = transform_to_json(
            transform
        )
        return image

    def transform_to_point_cloud(
        self, point_cloud_name: str, *, subcollection_name: str = "obsl"
    ) -> CoordinateTransform:
        """Returns the coordinate transform from the scene to a point cloud
        registered in the scene.


        TODO: Finish doc string.
        """
        try:
            subcollection: Collection = self[subcollection_name]  # type: ignore
        except KeyError as ke:
            raise KeyError(
                f"No collection '{subcollection_name}' in this scene."
            ) from ke
        try:
            transform_json = subcollection.metadata[
                f"soma_scene_registry_{point_cloud_name}"
            ]
        except KeyError as ke:
            raise KeyError(
                f"No coordinate space registry for '{point_cloud_name}' in collection "
                f"'{subcollection_name}'."
            ) from ke
        return transform_from_json(transform_json)

    def transform_to_multiscale_image(
        self,
        image_name: str,
        *,
        subcollection_name: str = "img",
        level: Optional[Union[str, int]] = None,
    ) -> CoordinateTransform:
        """Returns the corodinate transform from the scene to an image collection
        registered in the scene.

        If the name or level of an image is provided, the transformation will be to
        the requested level instead of the reference level of the multiscale image.
        """
        try:
            subcollection: Collection = self[subcollection_name]  # type: ignore
        except KeyError as ke:
            raise KeyError(
                f"No collection '{subcollection_name}' in this scene."
            ) from ke
        try:
            transform_json = subcollection.metadata[f"soma_scene_registry_{image_name}"]
        except KeyError:
            raise KeyError(
                f"No coordinate space registry for '{image_name}' in collection "
                f"'{subcollection_name}'"
            )
        base_transform = transform_from_json(transform_json)
        if level is None:
            return base_transform
        try:
            image: MultiscaleImage = subcollection[image_name]
        except KeyError as ke:
            raise KeyError(
                f"No MultiscaleImage named '{image_name}' in '{subcollection_name}'."
            ) from ke
        if isinstance(level, str):
            raise NotImplementedError(
                "Support for querying image level by name is not yet implemented."
            )
        level_transform = image.get_transformation_to_level(level)
        return base_transform * level_transform  # type: ignore[no-any-return]
