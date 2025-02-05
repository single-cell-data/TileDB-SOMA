# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.
"""
Implementation of a SOMA Scene
"""

from __future__ import annotations

import warnings
from typing import Any, List, Sequence, Tuple, Type, TypeVar, Union

import somacore
from somacore import (
    CoordinateSpace,
    CoordinateTransform,
    options,
)
from typing_extensions import Self

from . import _funcs, _tdb_handles
from . import pytiledbsoma as clib
from ._collection import CollectionBase
from ._constants import (
    SOMA_COORDINATE_SPACE_METADATA_KEY,
    SPATIAL_DISCLAIMER,
)
from ._exception import SOMAError, map_exception_for_create
from ._geometry_dataframe import GeometryDataFrame
from ._multiscale_image import MultiscaleImage
from ._point_cloud_dataframe import PointCloudDataFrame
from ._soma_object import AnySOMAObject
from ._spatial_util import (
    coordinate_space_from_json,
    coordinate_space_to_json,
    transform_from_json,
    transform_to_json,
)
from ._types import OpenTimestamp
from .options import SOMATileDBContext
from .options._soma_tiledb_context import _validate_soma_tiledb_context

_spatial_element = Union[GeometryDataFrame, MultiscaleImage, PointCloudDataFrame]

_SE = TypeVar("_SE", bound=_spatial_element)


class Scene(  # type: ignore[misc]   # __eq__ false positive
    CollectionBase[AnySOMAObject],
    somacore.Scene[
        MultiscaleImage, PointCloudDataFrame, GeometryDataFrame, AnySOMAObject
    ],
):
    """A collection subtype representing spatial assets that can all be stored
    on a single coordinate space.

    Lifecycle:
        Experimental.
    """

    __slots__ = ("_coord_space",)
    _wrapper_type = _tdb_handles.SceneWrapper

    _subclass_constrained_soma_types = {
        "img": ("SOMACollection",),
        "obsl": ("SOMACollection",),
        "varl": ("SOMACollection",),
    }

    @classmethod
    def create(
        cls,
        uri: str,
        *,
        coordinate_space: Sequence[str] | CoordinateSpace | None = None,
        platform_config: options.PlatformConfig | None = None,
        context: SOMATileDBContext | None = None,
        tiledb_timestamp: OpenTimestamp | None = None,
    ) -> Self:
        """Creates a new scene at the given URI.

        Args:
            uri:
                The location to create this SOMA scene at.
            coordinate_space:
                Optional coordinate space or the axis names for the coordinate space
                the scene is defined on. If ``None`` no coordinate space is set.
                Defaults to ``None``.
            platform_config:
                Platform-specific options used to create this scene. This may be
                provided as settings in a dictionary, with options located in the
                ``{'tiledb': {'create': ...}}`` key, or as a
                :class:`~tiledbsoma.TileDBCreateOptions` object.
            context:
                If provided, the :class:`SOMATileDBContext` to use when creating and
                opening this scene
            tiledb_timestamp:
                If specified, overrides the default timestamp used to open this object.
                If unset, uses the timestamp provided by the context.

        Returns:
            The newly created scene, opened for writing.

        Lifecycle:
            Experimental.
        """
        warnings.warn(SPATIAL_DISCLAIMER)

        context = _validate_soma_tiledb_context(context)

        if coordinate_space is None:
            axis_names = None
            axis_units = None
        elif isinstance(coordinate_space, CoordinateSpace):
            axis_names = tuple(axis.name for axis in coordinate_space)
            axis_units = tuple(axis.unit for axis in coordinate_space)
        else:
            axis_names = tuple(coordinate_space)
            axis_units = None

        try:
            timestamp_ms = context._open_timestamp_ms(tiledb_timestamp)
            clib.SOMAScene.create(
                ctx=context.native_context,
                uri=uri,
                axis_names=axis_names,
                axis_units=axis_units,
                timestamp=(0, timestamp_ms),
            )
            return cls(
                cls._wrapper_type.open(uri, "w", context, tiledb_timestamp),
                _dont_call_this_use_create_or_open_instead="tiledbsoma-internal-code",
            )
        except SOMAError as e:
            raise map_exception_for_create(e, uri) from None

    def __init__(
        self,
        handle: _tdb_handles.SOMAGroupWrapper[Any],
        **kwargs: Any,
    ):
        super().__init__(handle, **kwargs)
        coord_space = self.metadata.get(SOMA_COORDINATE_SPACE_METADATA_KEY)
        if coord_space is None:
            self._coord_space: CoordinateSpace | None = None
        else:
            self._coord_space = coordinate_space_from_json(coord_space)

    def _open_subcollection(
        self, subcollection: Union[str, Sequence[str]]
    ) -> CollectionBase[AnySOMAObject]:
        if len(subcollection) == 0:
            raise ValueError("Invalid subcollection: value cannot be empty.")
        if isinstance(subcollection, str):
            subcollection = (subcollection,)
        else:
            subcollection = tuple(subcollection)
        coll: CollectionBase[AnySOMAObject] = self
        # Keep track of collection hierarchy for informative error reporting
        parent_name: List[str] = []
        for name in subcollection:
            try:
                coll = coll[name]  # type: ignore[assignment]
            except KeyError as ke:
                raise KeyError(
                    f"Unable to open collection '{name}' in {parent_name}."
                ) from ke
            parent_name.append(name)
        return coll

    def _set_transform_to_element(
        self,
        kind: Type[_SE],
        *,
        key: str,
        transform: CoordinateTransform,
        subcollection: Union[str, Sequence[str]],
        coordinate_space: CoordinateSpace | None,
    ) -> _SE:
        # Check the transform is compatible with the coordinate spaces of the scene
        # and the new element coordinate space (if provided).
        if self.coordinate_space is None:
            raise SOMAError(
                "The scene coordinate space must be set before setting a transform."
            )
        if transform.input_axes != self.coordinate_space.axis_names:
            raise ValueError(
                f"The name of the transform input axes, {transform.input_axes}, do "
                f"not match the name of the axes, {self.coordinate_space.axis_names}, "
                f"in the scene coordinate space."
            )
        if (
            coordinate_space is not None
            and transform.output_axes != coordinate_space.axis_names
        ):
            raise ValueError(
                f"The name of the transform output axes, {transform.output_axes}, do "
                f"not match the name of the axes, {coordinate_space.axis_names}, ."
                f" in the provided coordinate space."
            )

        # Check asset exists in the specified location.
        coll = self._open_subcollection(subcollection)
        try:
            elem = coll[key]
        except KeyError as ke:
            raise KeyError(f"No element named '{key}' in '{subcollection}'.") from ke
        if not isinstance(elem, kind):
            raise TypeError(
                f"'{key}' in '{subcollection}' is a {type(elem).__name__} not a {kind.__name__}."
            )

        # Either set the new coordinate space or check the axes of the current
        # coordinate space the element is defined on.
        if coordinate_space is None:
            elem_axis_names: Tuple[str, ...] = elem.coordinate_space.axis_names
            if elem_axis_names != transform.output_axes:
                raise ValueError(
                    f"The name of transform output axes, {transform.output_axes}, do "
                    f"not match the name of the axes in the multiscale image coordinate"
                    f" space, {elem_axis_names}."
                )
        else:
            elem.coordinate_space = coordinate_space

        # Set the transform metadata and return the multisclae image.
        coll.metadata[f"soma_scene_registry_{key}"] = transform_to_json(transform)
        return elem

    @property
    def coordinate_space(self) -> CoordinateSpace | None:
        """Coordinate system for this scene.

        Lifecycle:
            Experimental.
        """
        return self._coord_space

    @coordinate_space.setter
    def coordinate_space(self, value: CoordinateSpace) -> None:
        if not isinstance(value, CoordinateSpace):
            raise TypeError(f"Invalid type {type(value).__name__}.")
        self.metadata[SOMA_COORDINATE_SPACE_METADATA_KEY] = coordinate_space_to_json(
            value
        )
        self._coord_space = value

    @_funcs.forwards_kwargs_to(
        GeometryDataFrame.create, exclude=("context", "tiledb_timestamp")
    )
    def add_new_geometry_dataframe(
        self,
        key: str,
        subcollection: Union[str, Sequence[str]],
        *,
        transform: CoordinateTransform | None,
        uri: str | None = None,
        **kwargs: Any,
    ) -> GeometryDataFrame:
        """Adds a ``GeometryDataFrame`` to the scene and sets a coordinate transform
        between the scene and the dataframe.

        If the subcollection the geometry dataframe will be created inside of is more
        than one layer deep, the input should be provided as a sequence of names. For
        example, to add a new geometry dataframe named  "transcripts" in the "var/RNA"
        collection::

            scene.add_new_geometry_dataframe(
                'transcripts', subcollection=['var', 'RNA'], **kwargs
            )

        See :meth:`add_new_collection` for details about child URIs.

        Args:
            key: The name of the geometry dataframe.
            subcollection: The name, or sequence of names, of the subcollection the
                dataframe is stored in.
            transform: The coordinate transformation from the scene to the dataframe.
            uri: If provided, overrides the default URI what would be used to create
                this object. This may be aboslution or relative.
            kwargs: Additional keyword arugments as specified in
                :meth:`spatial.GeometryDataFrame.create`.

        Returns:
            The newly create ``GeometryDataFrame``, opened for writing.

        Lifecycle:
            Experimental.
        """
        raise NotImplementedError()

    @_funcs.forwards_kwargs_to(
        MultiscaleImage.create, exclude=("context", "tiledb_timestamp")
    )
    def add_new_multiscale_image(
        self,
        key: str,
        subcollection: Union[str, Sequence[str]],
        *,
        transform: CoordinateTransform | None,
        uri: str | None = None,
        coordinate_space: Union[Sequence[str], CoordinateSpace] = ("x", "y"),
        **kwargs: Any,
    ) -> MultiscaleImage:
        """Adds a ``MultiscaleImage`` to the scene and sets a coordinate transform
        between the scene and the dataframe.

        See :meth:`add_new_collection` for details about child URIs.

        Args:
            key: The name of the multiscale image.
            subcollection: The name, or sequence of names, of the subcollection the
                dataframe is stored in.
            transform: The coordinate transformation from the scene to the dataframe.
            uri: If provided, overrides the default URI what would be used to create
                this object. This may be aboslution or relative.
            kwargs: Additional keyword arugments as specified in
                :meth:`spatial.MultiscaleImage.create`.

        Returns:
            The newly create ``MultiscaleImage``, opened for writing.

        Lifecycle:
            Experimental.
        """
        if transform is not None:
            # Get and check the scene coordinate space axis names.
            if self.coordinate_space is None:
                raise SOMAError(
                    "The scene coordinate space must be set before setting a transform."
                )
            if transform.input_axes != self.coordinate_space.axis_names:
                raise ValueError(
                    f"The name of the transform input axes, {transform.input_axes}, "
                    f"do not match the name of the axes, "
                    f"{self.coordinate_space.axis_names}, in the scene coordinate "
                    f"space."
                )

            if transform.input_axes != self.coordinate_space.axis_names:
                raise ValueError(
                    f"The name of the transform input axes, {transform.input_axes}, "
                    f"do not match the name of the axes, "
                    f"{self.coordinate_space.axis_names}, in the scene coordinate "
                    f"space."
                )

            # Get multisclae image coordinate space and check.
            elem_axis_names = (
                coordinate_space.axis_names
                if isinstance(coordinate_space, CoordinateSpace)
                else tuple(coordinate_space)
            )
            if transform.output_axes != elem_axis_names:
                raise ValueError(
                    f"The name of the transform output axes, {transform.output_axes}, "
                    f"do not match the name of the axes, {elem_axis_names}, of the "
                    f"coordinate space the multiscale image is defined on."
                )

        # Open the subcollection and add the new multiscale image.
        coll = self._open_subcollection(subcollection)
        image = coll._add_new_element(
            key,
            MultiscaleImage,
            lambda create_uri: MultiscaleImage.create(
                create_uri,
                context=self.context,
                tiledb_timestamp=self.tiledb_timestamp_ms,
                coordinate_space=coordinate_space,
                **kwargs,
            ),
            uri,
        )

        # Store the metadata for the transform.
        if transform is not None:
            coll.metadata[f"soma_scene_registry_{key}"] = transform_to_json(transform)

        # Return the multiscale image.
        return image

    @_funcs.forwards_kwargs_to(
        PointCloudDataFrame.create, exclude=("context", "tiledb_timestamp")
    )
    def add_new_point_cloud_dataframe(
        self,
        key: str,
        subcollection: Union[str, Sequence[str]],
        *,
        transform: CoordinateTransform | None,
        uri: str | None = None,
        coordinate_space: Union[Sequence[str], CoordinateSpace] = ("x", "y"),
        **kwargs: Any,
    ) -> PointCloudDataFrame:
        """Adds a point cloud to the scene and sets a coordinate transform
        between the scene and the dataframe.

        If the subcollection the point cloud dataframe will be added to is more than
        one layer deep, the input should be provided as a sequence of names. For
        example, to add a new point cloud dataframe named  "transcripts" to the
        "var/RNA" collection::

            scene.add_new_point_cloud_dataframe(
                'transcripts', subcollection=['var', 'RNA'], **kwargs
            )


        See :meth:`add_new_collection` for details about child URIs.

        Args:
            key: The name of the point cloud dataframe.
            subcollection: The name, or sequence of names, of the subcollection the
                dataframe is stored in.
            transform: The coordinate transformation from the scene to the dataframe.
            uri: If provided, overrides the default URI what would be used to create
                this object. This may be aboslution or relative.
            kwargs: Additional keyword arugments as specified in
                :meth:`spatial.PointCloudDataFrame.create`.

        Returns:
            The newly created ``PointCloudDataFrame``, opened for writing.

        Lifecycle:
            Experimental.
        """
        # If the transform is set, check it is consistent with the coordinate spaces.
        if transform is not None:
            # Get Scene coordinate space and check the axis names.
            if self.coordinate_space is None:
                raise SOMAError(
                    "The scene coordinate space must be set before setting a transform."
                )
            if transform.input_axes != self.coordinate_space.axis_names:
                raise ValueError(
                    f"The name of the transform input axes, {transform.input_axes}, "
                    f"do not match the name of the axes, "
                    f"{self.coordinate_space.axis_names}, in the scene coordinate "
                    f"space."
                )

            # Get point cloud coordinate space and check
            elem_axis_names = (
                coordinate_space.axis_names
                if isinstance(coordinate_space, CoordinateSpace)
                else tuple(coordinate_space)
            )
            if transform.output_axes != elem_axis_names:
                raise ValueError(
                    f"The name of the transform output axes, {transform.output_axes}, "
                    f"do not match the name of the axes, {elem_axis_names}, of the "
                    f"coordinate space the point cloud is defined on."
                )

        # Open the collection and add the new point cloud.
        coll = self._open_subcollection(subcollection)
        point_cloud = coll._add_new_element(
            key,
            PointCloudDataFrame,
            lambda create_uri: PointCloudDataFrame.create(
                create_uri,
                context=self.context,
                tiledb_timestamp=self.tiledb_timestamp_ms,
                coordinate_space=coordinate_space,
                **kwargs,
            ),
            uri,
        )

        # Store the metadata for the transform.
        if transform is not None:
            coll.metadata[f"soma_scene_registry_{key}"] = transform_to_json(transform)

        # Return the point cloud.
        return point_cloud

    def set_transform_to_geometry_dataframe(
        self,
        key: str,
        subcollection: Union[str, Sequence[str]] = "obsl",
        *,
        transform: CoordinateTransform,
        coordinate_space: CoordinateSpace | None = None,
    ) -> GeometryDataFrame:
        """Adds the coordinate transform for the scene coordinate space to
        a geometry dataframe stored in the scene.

        If the subcollection the geometry dataframe is inside of is more than one
        layer deep, the input should be provided as a sequence of names. For example,
        to set a transformation for geometry dataframe named  "transcripts" in the
        "var/RNA" collection::

            scene.set_transform_to_geometry_dataframe(
                'transcripts', transform, subcollection=['var', 'RNA'],
            )

        Args:
            key: The name of the geometry dataframe.
            subcollection: The name, or sequence of names, of the subcollection the
                dataframe is stored in. Defaults to ``'obsl'``.
            transform: The coordinate transformation from the scene to the dataframe.
            coordinate_space: Optional coordinate space for the dataframe. This will
                replace the existing coordinate space of the dataframe.

        Returns:
            The geometry dataframe, opened for writing.

        Lifecycle: experimental
        """
        raise NotImplementedError()

    def set_transform_to_multiscale_image(
        self,
        key: str,
        subcollection: Union[str, Sequence[str]] = "img",
        *,
        transform: CoordinateTransform,
        coordinate_space: CoordinateSpace | None = None,
    ) -> MultiscaleImage:
        """Adds the coordinate transform for the scene coordinate space to
        a multiscale image stored in the scene.

        The transform to the multiscale image must be to the coordinate space
        defined on the reference level for the image. In most cases, this will be
        the level ``0`` image.

        Args:
            key: The name of the multiscale image.
            subcollection: The name, or sequence of names, of the subcollection the
                image is stored in. Defaults to ``'img'``.
            transform: The coordinate transformation from the scene to the reference
                level of the multiscale image.
            coordinate_space: Optional coordinate space for the image. This will
                replace the existing coordinate space of the multiscale image.

        Returns:
            The multiscale image, opened for writing.

        Lifecycle: experimental
        """
        return self._set_transform_to_element(
            MultiscaleImage,
            key=key,
            transform=transform,
            subcollection=subcollection,
            coordinate_space=coordinate_space,
        )

    def set_transform_to_point_cloud_dataframe(
        self,
        key: str,
        subcollection: Union[str, Sequence[str]] = "obsl",
        *,
        transform: CoordinateTransform,
        coordinate_space: CoordinateSpace | None = None,
    ) -> PointCloudDataFrame:
        """Adds the coordinate transform for the scene coordinate space to
        a point cloud dataframe stored in the scene.

        If the subcollection the point cloud dataframe is inside of is more than one
        layer deep, the input should be provided as a sequence of names. For example,
        to set a transform for  a point named `transcripts` in the `var/RNA`
        collection:

            scene.set_transform_to_point_cloud_dataframe(
                'transcripts', transform, subcollection=['var', 'RNA'],
            )

        Args:
            key: The name of the point cloud dataframe.
            transform: The coordinate transformation from the scene to the dataframe.
            subcollection: The name, or sequence of names, of the subcollection the
                dataframe is stored in. Defaults to ``'obsl'``.
            coordinate_space: Optional coordinate space for the dataframe. This will
                replace the existing coordinate space of the dataframe. Defaults to
                ``None``.

        Returns:
            The point cloud dataframe, opened for writing.

        Lifecycle: experimental
        """
        return self._set_transform_to_element(
            PointCloudDataFrame,
            key=key,
            transform=transform,
            subcollection=subcollection,
            coordinate_space=coordinate_space,
        )

    def get_transform_from_geometry_dataframe(
        self, key: str, subcollection: Union[str, Sequence[str]] = "obsl"
    ) -> CoordinateTransform:
        """Returns the coordinate transformation from the requested geometry dataframe
        to the scene.

        Args:
            key: The name of the geometry dataframe.
            subcollection: The name, or sequence of names, of the subcollection the
                dataframe is stored in. Defaults to ``'obsl'``.

        Returns:
            Coordinate transform from the dataframe to the scene.

        Lifecycle: experimental
        """
        transform = self.get_transform_to_geometry_dataframe(key, subcollection)
        return transform.inverse_transform()

    def get_transform_from_multiscale_image(
        self,
        key: str,
        subcollection: str = "img",
        *,
        level: str | int | None = None,
    ) -> CoordinateTransform:
        """Returns the coordinate transformation from the requested multiscale image to
        the scene.

        Args:
            key: The name of the multiscale image.
            subcollection: The name, or sequence of names, of the subcollection the
                dataframe is stored in. Defaults to ``'img'``.
            level: The level of the image to get the transformation from.
                Defaults to ``None`` -- the transformation will be to the reference
                level.

        Returns:
            Coordinate transform from the multiscale image to the scene.

        Lifecycle: experimental
        """
        if level is None:
            transform = self.get_transform_to_multiscale_image(key, subcollection)
            return transform.inverse_transform()
        coll = self._open_subcollection(subcollection)
        try:
            transform_json = coll.metadata[f"soma_scene_registry_{key}"]
        except KeyError:
            raise KeyError(
                f"No coordinate space registry for '{key}' in collection "
                f"'{subcollection}'"
            )
        base_transform = transform_from_json(transform_json)
        try:
            image: MultiscaleImage = coll[key]  # type: ignore[assignment]
        except KeyError as ke:
            raise KeyError(
                f"No MultiscaleImage named '{key}' in '{subcollection}'."
            ) from ke
        if not isinstance(image, MultiscaleImage):
            raise TypeError(
                f"Item at '{key}' in '{subcollection}' has an unexpected type "
                f"{type(image)!r}."
            )
        level_transform = image.get_transform_from_level(level)
        return base_transform.inverse_transform() @ level_transform

    def get_transform_from_point_cloud_dataframe(
        self, key: str, subcollection: str = "obsl"
    ) -> CoordinateTransform:
        """Returns the coordinate transformation from the requested point cloud
        dataframe to the scene.

        Args:
            key: The name of the point cloud dataframe.
            subcollection: The name, or sequence of names, of the subcollection the
                dataframe is stored in. Defaults to ``'obsl'``.

        Returns:
            Coordinate transform from the scene to the point cloud dataframe.

        Lifecycle: experimental
        """
        transform = self.get_transform_to_point_cloud_dataframe(
            key, subcollection=subcollection
        )
        return transform.inverse_transform()

    def get_transform_to_geometry_dataframe(
        self, key: str, subcollection: Union[str, Sequence[str]] = "obsl"
    ) -> CoordinateTransform:
        """Returns the coordinate transformation from the scene to a requested
        geometery dataframe.

        Args:
            key: The name of the geometry dataframe.
            subcollection: The name, or sequence of names, of the subcollection the
                dataframe is stored in. Defaults to ``'obsl'``.

        Returns:
            Coordinate transform from the scene to the requested dataframe.

        Lifecycle: experimental
        """
        coll = self._open_subcollection(subcollection)
        try:
            transform_json = coll.metadata[f"soma_scene_registry_{key}"]
        except KeyError as ke:
            raise KeyError(
                f"No coordinate space registry for '{key}' in collection "
                f"'{subcollection}'."
            ) from ke
        return transform_from_json(transform_json)

    def get_transform_to_multiscale_image(
        self,
        key: str,
        subcollection: str = "img",
        *,
        level: str | int | None = None,
    ) -> CoordinateTransform:
        """Returns the coordinate transformation from the scene to a requested
        multiscale image.

        Args:
            key: The name of the multiscale image.
            subcollection: The name, or sequence of names, of the subcollection the
                dataframe is stored in. Defaults to ``'img'``.
            level: The level of the image to get the transformation to.
                Defaults to ``None`` -- the transformation will be to the reference
                level.

        Returns:
            Coordinate transform from the scene to the requested multiscale image.

        Lifecycle: experimental
        """
        coll = self._open_subcollection(subcollection)
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
            image: MultiscaleImage = coll[key]  # type: ignore[assignment]
        except KeyError as ke:
            raise KeyError(
                f"No MultiscaleImage named '{key}' in '{subcollection}'."
            ) from ke
        if not isinstance(image, MultiscaleImage):
            raise TypeError(
                f"Item at '{key}' in '{subcollection}' has an unexpected type "
                f"{type(image)!r}."
            )
        level_transform = image.get_transform_to_level(level)
        return level_transform @ base_transform

    def get_transform_to_point_cloud_dataframe(
        self, key: str, subcollection: str = "obsl"
    ) -> CoordinateTransform:
        """Returns the coordinate transformation from the scene to a requested
        point cloud dataframe.

        Args:
            key: The name of the point cloud dataframe.
            subcollection: The name, or sequence of names, of the subcollection the
                dataframe is stored in. Defaults to ``'obsl'``.

        Returns:
            Coordinate transform from the scene to the point cloud dataframe.

        Lifecycle: experimental
        """
        coll = self._open_subcollection(subcollection)
        try:
            transform_json = coll.metadata[f"soma_scene_registry_{key}"]
        except KeyError as ke:
            raise KeyError(
                f"No coordinate space registry for '{key}' in collection "
                f"'{subcollection}'."
            ) from ke
        return transform_from_json(transform_json)
