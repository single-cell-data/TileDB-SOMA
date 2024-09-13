# Copyright (c) 2024 TileDB, Inc.
#
# Licensed under the MIT License.

"""Implementation of a SOMA Scene."""

from typing import Any, Optional, Sequence, Union

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
    scene.Scene[MultiscaleImage, PointCloud, GeometryDataFrame, AnySOMAObject],
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

    def register_geometry_dataframe(
        self,
        key: str,
        transform: CoordinateTransform,
        *,
        subcollection: Union[str, Sequence[str]] = "obsl",
        coordinate_space: Optional[CoordinateSpace] = None,
    ) -> GeometryDataFrame:
        raise NotImplementedError(
            "Support for registering geometry data frames is not yet implemented."
        )

    def register_multiscale_image(
        self,
        key: str,
        transform: CoordinateTransform,
        *,
        subcollection: Union[str, Sequence[str]] = "img",
        coordinate_space: Optional[CoordinateSpace] = None,
    ) -> MultiscaleImage:
        """TODO Add docstring"""
        if not isinstance(subcollection, str):
            raise NotImplementedError()

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
            coll: Collection = self[subcollection]  # type: ignore
        except KeyError as ke:
            raise KeyError(f"No collection '{subcollection}' in this scene.") from ke
        try:
            image: MultiscaleImage = coll[key]
        except KeyError as ke:
            raise KeyError(
                f"No multiscale image named '{key}' in '{subcollection}'."
            ) from ke
        if not isinstance(image, MultiscaleImage):
            raise TypeError(f"'{key}' in '{subcollection}' is not an MultiscaleImage.")

        image.coordinate_space = coordinate_space
        coll.metadata[f"soma_scene_registry_{key}"] = transform_to_json(transform)
        return image

    def register_point_cloud(
        self,
        key: str,
        transform: CoordinateTransform,
        *,
        subcollection: Union[str, Sequence[str]] = "obsl",
        coordinate_space: Optional[CoordinateSpace] = None,
    ) -> PointCloud:
        if not isinstance(subcollection, str):
            raise NotImplementedError()
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
            coll: Collection = self[subcollection]  # type: ignore
        except KeyError as ke:
            raise KeyError(f"No collection '{subcollection}' in this scene.") from ke
        try:
            point_cloud: PointCloud = coll[key]
        except KeyError as ke:
            raise KeyError(f"No PointCloud named '{key}' in '{coll}'.") from ke
        if not isinstance(point_cloud, PointCloud):
            raise TypeError(f"'{key}' in '{subcollection}' is not an PointCloud.")

        point_cloud.coordinate_space = coordinate_space
        coll.metadata[f"soma_scene_registry_{key}"] = transform_to_json(transform)
        return point_cloud

    def get_transformation_to_geometry_dataframe(
        self, key: str, *, subcollection: Union[str, Sequence[str]] = "obsl"
    ) -> CoordinateTransform:
        """TODO: Add doc string."""
        if not isinstance(subcollection, str):
            raise NotImplementedError()
        try:
            coll: Collection = self[subcollection]  # type: ignore
        except KeyError as ke:
            raise KeyError(f"No collection '{subcollection}' in this scene.") from ke
        try:
            transform_json = coll.metadata[f"soma_scene_registry_{key}"]
        except KeyError as ke:
            raise KeyError(
                f"No coordinate space registry for '{key}' in collection "
                f"'{subcollection}'."
            ) from ke
        return transform_from_json(transform_json)

    def get_transformation_to_multiscale_image(
        self,
        key: str,
        *,
        subcollection: Union[str, Sequence[str]] = "img",
        level: Optional[Union[str, int]] = None,
    ) -> CoordinateTransform:
        """Returns the corodinate transform from the scene to an image collection
        registered in the scene.

        If the name or level of an image is provided, the transformation will be to
        the requested level instead of the reference level of the multiscale image.
        """
        if not isinstance(subcollection, str):
            raise NotImplementedError()
        try:
            coll: Collection = self[subcollection]  # type: ignore
        except KeyError as ke:
            raise KeyError(f"No collection '{subcollection}' in this scene.") from ke
        try:
            transform_json = coll.metadata[f"soma_scene_registry_{key}"]
        except KeyError:
            raise KeyError(
                f"No coordinate space registry for '{key}' in collection "
                f"'{subcollection}'"
            )
        base_transform = transform_from_json(transform_json)
        if level is None:
            return base_transform
        try:
            image: MultiscaleImage = coll[key]
        except KeyError as ke:
            raise KeyError(
                f"No MultiscaleImage named '{key}' in '{subcollection}'."
            ) from ke
        if isinstance(level, str):
            raise NotImplementedError(
                "Support for querying image level by name is not yet implemented."
            )
        level_transform = image.get_transformation_to_level(level)
        return level_transform @ base_transform

    def get_transformation_to_point_cloud(
        self, key: str, *, subcollection: Union[str, Sequence[str]] = "obsl"
    ) -> CoordinateTransform:
        """Returns the coordinate transform from the scene to a point cloud
        registered in the scene.


        TODO: Finish doc string.
        """
        if not isinstance(subcollection, str):
            raise NotImplementedError()
        try:
            coll: Collection = self[subcollection]  # type: ignore
        except KeyError as ke:
            raise KeyError(f"No collection '{subcollection}' in this scene.") from ke
        try:
            transform_json = coll.metadata[f"soma_scene_registry_{key}"]
        except KeyError as ke:
            raise KeyError(
                f"No coordinate space registry for '{key}' in collection "
                f"'{subcollection}'."
            ) from ke
        return transform_from_json(transform_json)
