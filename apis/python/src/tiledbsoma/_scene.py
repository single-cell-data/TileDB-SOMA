# Copyright (c) 2024 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2024 TileDB, Inc.
#
# Licensed under the MIT License.
"""
Implementation of a SOMA Scene
"""

from typing import Any, List, Optional, Sequence, Union

import somacore
from somacore import (
    CoordinateSpace,
    CoordinateTransform,
)

from . import _funcs, _tdb_handles
from ._collection import Collection, CollectionBase
from ._constants import SOMA_COORDINATE_SPACE_METADATA_KEY
from ._exception import SOMAError
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

    def _check_transform_to_asset(
        self,
        transform: CoordinateTransform,
        asset_coord_space: Optional[CoordinateSpace],
    ) -> None:
        """Raises an error if the scene coordinate space is not set or the
        scene axis names do not match the input transform axis names.
        """
        if self.coordinate_space is None:
            raise SOMAError(
                "The scene coordinate space must be set before setting a transform."
            )
        if transform.input_axes != self.coordinate_space.axis_names:
            raise ValueError(
                f"The name of the transform input axes, {transform.input_axes}, do "
                f"not match the name of the axes in the scene coordinate space, "
                f"{self.coordinate_space.axis_names}."
            )
        if (
            asset_coord_space is not None
            and transform.output_axes != asset_coord_space.axis_names
        ):
            raise ValueError(
                f"The name of the transform output axes, {transform.output_axes}, do "
                f"not match the name of the axes in the provided coordinate space, "
                f"{asset_coord_space.axis_names}."
            )

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

    @property
    def coordinate_space(self) -> Optional[CoordinateSpace]:
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
        transform: Optional[CoordinateTransform],
        uri: Optional[str] = None,
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
                dataframe is stored in. Defaults to ``'obsl'``.
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
        transform: Optional[CoordinateTransform],
        uri: Optional[str] = None,
        **kwargs: Any,
    ) -> MultiscaleImage:
        """Adds a ``MultiscaleImage`` to the scene and sets a coordinate transform
        between the scene and the dataframe.

        See :meth:`add_new_collection` for details about child URIs.

        Args:
            key: The name of the multiscale image.
            subcollection: The name, or sequence of names, of the subcollection the
                dataframe is stored in. Defaults to ``'obsl'``.
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
            raise NotImplementedError()

        coll = self._open_subcollection(subcollection)
        return coll._add_new_element(
            key,
            MultiscaleImage,
            lambda create_uri: MultiscaleImage.create(
                create_uri,
                context=self.context,
                tiledb_timestamp=self.tiledb_timestamp_ms,
                **kwargs,
            ),
            uri,
        )

    @_funcs.forwards_kwargs_to(
        PointCloudDataFrame.create, exclude=("context", "tiledb_timestamp")
    )
    def add_new_point_cloud_dataframe(
        self,
        key: str,
        subcollection: Union[str, Sequence[str]],
        *,
        transform: Optional[CoordinateTransform],
        uri: Optional[str] = None,
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
                dataframe is stored in. Defaults to ``'obsl'``.
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
        if transform is not None:
            raise NotImplementedError()
        coll = self._open_subcollection(subcollection)
        return coll._add_new_element(
            key,
            PointCloudDataFrame,
            lambda create_uri: PointCloudDataFrame.create(
                create_uri,
                context=self.context,
                tiledb_timestamp=self.tiledb_timestamp_ms,
                **kwargs,
            ),
            uri,
        )

    def set_transform_to_geometry_dataframe(
        self,
        key: str,
        transform: CoordinateTransform,
        *,
        subcollection: Union[str, Sequence[str]] = "obsl",
        coordinate_space: Optional[CoordinateSpace] = None,
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
            transform: The coordinate transformation from the scene to the dataframe.
            subcollection: The name, or sequence of names, of the subcollection the
                dataframe is stored in. Defaults to ``'obsl'``.
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
        transform: CoordinateTransform,
        *,
        subcollection: Union[str, Sequence[str]] = "img",
        coordinate_space: Optional[CoordinateSpace] = None,
    ) -> MultiscaleImage:
        """Adds the coordinate transform for the scene coordinate space to
        a multiscale image stored in the scene.

        The transform to the multiscale image must be to the coordinate space
        defined on the reference level for the image. In most cases, this will be
        the level ``0`` image.

        Args:
            key: The name of the multiscale image.
            transform: The coordinate transformation from the scene to the reference
                level of the multiscale image.
            subcollection: The name, or sequence of names, of the subcollection the
                image is stored in. Defaults to ``'img'``.
            coordinate_space: Optional coordinate space for the image. This will
                replace the existing coordinate space of the multiscale image.

        Returns:
            The multiscale image, opened for writing.

        Lifecycle: experimental
        """
        # Check the transform is compatible with the oordinate spaces.
        self._check_transform_to_asset(transform, coordinate_space)
        assert (
            self.coordinate_space is not None
        )  # Assert for typing - verified in the above method.

        # Check asset exists in the specified location.
        coll = self._open_subcollection(subcollection)
        try:
            image: MultiscaleImage = coll[key]  # type: ignore[assignment]
        except KeyError as ke:
            raise KeyError(
                f"No multiscale image named '{key}' in '{subcollection}'."
            ) from ke
        if not isinstance(image, MultiscaleImage):
            raise TypeError(f"'{key}' in '{subcollection}' is not an MultiscaleImage.")

        # Either set the new coordinate space or check the axes of the current
        # coordinate space the multiscale image is defined on.
        if coordinate_space is None:
            if image.coordinate_space.axis_names != transform.output_axes:
                raise ValueError(
                    f"The name of transform output axes, {transform.output_axes}, do"
                    f"not match the name of the axes in the multiscale image coordinate"
                    f" space, {image.coordinate_space.axis_names}."
                )
        else:
            image.coordinate_space = coordinate_space

        # Set the transform metadata and return the multisclae image.
        coll.metadata[f"soma_scene_registry_{key}"] = transform_to_json(transform)
        return image

    def set_transform_to_point_cloud_dataframe(
        self,
        key: str,
        transform: CoordinateTransform,
        *,
        subcollection: Union[str, Sequence[str]] = "obsl",
        coordinate_space: Optional[CoordinateSpace] = None,
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
        # Check the transform is compatible with the Scene coordinate spaces.
        self._check_transform_to_asset(transform, coordinate_space)
        assert (
            self.coordinate_space is not None
        )  # Assert for typing - verified in the above method.

        # Check asset exists in the specified location.
        try:
            coll: Collection = self[subcollection]  # type: ignore
        except KeyError as ke:
            raise KeyError(f"No collection '{subcollection}' in this scene.") from ke
        try:
            point_cloud: PointCloudDataFrame = coll[key]
        except KeyError as ke:
            raise KeyError(f"No PointCloudDataFrame named '{key}' in '{coll}'.") from ke
        if not isinstance(point_cloud, PointCloudDataFrame):
            raise TypeError(
                f"'{key}' in '{subcollection}' is not an PointCloudDataFrame."
            )

        # Either set the new coordinate space or check the axes of the current point
        # cloud coordinate space.
        if coordinate_space is None:
            if point_cloud.coordinate_space.axis_names != transform.output_axes:
                raise ValueError(
                    f"The name of transform output axes, {transform.output_axes}, do"
                    f"not match the name of the axes in the point cloud coordinate "
                    f"space, {point_cloud.coordinate_space.axis_names}."
                )
        else:
            point_cloud.coordinate_space = coordinate_space

        # Set the transform metadata and return the point cloud.
        coll.metadata[f"soma_scene_registry_{key}"] = transform_to_json(transform)
        return point_cloud

    def get_transform_from_geometry_dataframe(
        self, key: str, *, subcollection: Union[str, Sequence[str]] = "obsl"
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
        raise NotImplementedError()

    def get_transform_from_multiscale_image(
        self,
        key: str,
        *,
        subcollection: str = "img",
        level: Optional[Union[str, int]] = None,
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
        raise NotImplementedError()

    def get_transform_from_point_cloud_dataframe(
        self, key: str, *, subcollection: str = "obsl"
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
        raise NotImplementedError()

    def get_transform_to_geometry_dataframe(
        self, key: str, *, subcollection: Union[str, Sequence[str]] = "obsl"
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
        *,
        subcollection: str = "img",
        level: Optional[Union[str, int]] = None,
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
        self, key: str, *, subcollection: str = "obsl"
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
