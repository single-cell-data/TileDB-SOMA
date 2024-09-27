# Copyright (c) 2024 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2024 TileDB, Inc.
#
# Licensed under the MIT License.
"""
Implementation of a SOMA Scene
"""

from typing import Any, Optional, Sequence, Tuple, Union

import pyarrow as pa
import somacore
from somacore import (
    CoordinateSpace,
    CoordinateTransform,
    options,
)

from ._collection import CollectionBase
from ._constants import SOMA_GEOMETRY, SOMA_JOINID
from ._geometry_dataframe import GeometryDataFrame
from ._multiscale_image import MultiscaleImage
from ._point_cloud import PointCloud
from ._soma_object import AnySOMAObject


class Scene(  # type: ignore[misc]   # __eq__ false positive
    CollectionBase[AnySOMAObject],
    somacore.Scene[MultiscaleImage, PointCloud, GeometryDataFrame, AnySOMAObject],
):
    """A collection subtype representing spatial assets that can all be stored
    on a single coordinate space.

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
    def coordinate_space(self) -> Optional[CoordinateSpace]:
        """Coordinate system for this scene.

        Lifecycle:
            Experimental.
        """
        raise NotImplementedError()

    @coordinate_space.setter
    def coordinate_space(self, value: CoordinateSpace) -> None:
        raise NotImplementedError()

    def add_geometry_dataframe(
        self,
        key: str,
        subcollection: Union[str, Sequence[str]],
        transform: Optional[CoordinateTransform],
        *,
        uri: str,
        schema: pa.Schema,
        index_column_names: Sequence[str] = (SOMA_JOINID, SOMA_GEOMETRY),
        axis_names: Sequence[str] = ("x", "y"),
        domain: Optional[Sequence[Optional[Tuple[Any, Any]]]] = None,
        platform_config: Optional[options.PlatformConfig] = None,
        context: Optional[Any] = None,
    ) -> GeometryDataFrame:
        """Adds a ``GeometryDataFrame`` to the scene and sets a coordinate transform
        between the scene and the dataframe.

        If the subcollection the geometry dataframe is inside of is more than one
        layer deep, the input should be provided as a sequence of names. For example,
        to set the transformation to a geometry dataframe named  "transcripts" in
        the "var/RNA" collection::

            scene.add_geometry_dataframe(
                'cell_boundaries', subcollection=['var', 'RNA'], **kwargs
            )

        Args:
            key: The name of the geometry dataframe.
            transform: The coordinate transformation from the scene to the dataframe.
            subcollection: The name, or sequence of names, of the subcollection the
                dataframe is stored in. Defaults to ``'obsl'``.

        Returns:
            The newly create ``GeometryDataFrame``, opened for writing.

        Lifecycle: experimental
        """
        raise NotImplementedError()

    def add_multiscale_image(
        self,
        key: str,
        subcollection: Union[str, Sequence[str]],
        transform: Optional[CoordinateTransform],
        *,
        uri: str,
        type: pa.DataType,
        image_type: str = "CYX",
        reference_level_shape: Sequence[int],
        axis_names: Sequence[str] = ("c", "x", "y"),
    ) -> MultiscaleImage:
        """Adds a ``MultiscaleImage`` to the scene and sets a coordinate transform
        between the scene and the dataframe.

        Parameters are as in :meth:`spatial.PointCloud.create`.
        See :meth:`add_new_collection` for details about child URIs.

        Args:
            key: The name of the geometry dataframe.
            transform: The coordinate transformation from the scene to the dataframe.
            subcollection: The name, or sequence of names, of the subcollection the
                dataframe is stored in. Defaults to ``'obsl'``.

        Returns:
            The newly create ``MultiscaleImage``, opened for writing.

        Lifecycle: experimental
        """
        raise NotImplementedError()

    def add_new_point_cloud(
        self,
        key: str,
        subcollection: Union[str, Sequence[str]],
        transform: Optional[CoordinateTransform],
        *,
        uri: Optional[str] = None,
        schema: pa.Schema,
        index_column_names: Sequence[str] = (SOMA_JOINID,),
        axis_names: Sequence[str] = ("x", "y"),
        domain: Optional[Sequence[Optional[Tuple[Any, Any]]]] = None,
        platform_config: Optional[options.PlatformConfig] = None,
    ) -> PointCloud:
        """Adds a point cloud to the scene and sets a coordinate transform
        between the scene and the dataframe.

        Parameters are as in :meth:`spatial.PointCloud.create`.
        See :meth:`add_new_collection` for details about child URIs.

        Args:
            key: The name of the geometry dataframe.
            transform: The coordinate transformation from the scene to the dataframe.
            subcollection: The name, or sequence of names, of the subcollection the
                dataframe is stored in. Defaults to ``'obsl'``.

        Returns:
            The newly created ``PointCloud``, opened for writing.

        Lifecycle: experimental
        """
        raise NotImplementedError()

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

            scene.set_transfrom_for_geometry_dataframe(
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
        raise NotImplementedError()

    def set_transform_to_point_cloud(
        self,
        key: str,
        transform: CoordinateTransform,
        *,
        subcollection: Union[str, Sequence[str]] = "obsl",
        coordinate_space: Optional[CoordinateSpace] = None,
    ) -> PointCloud:
        """Adds the coordinate transform for the scene coordinate space to
        a point cloud stored in the scene.

        If the subcollection the point cloud is inside of is more than one
        layer deep, the input should be provided as a sequence of names. For example,
        to set a transform for  a point named `transcripts` in the `var/RNA`
        collection::

            scene.set_transformation_for_point_cloud(
                'transcripts', transform, subcollection=['var', 'RNA'],
            )

        Args:
            key: The name of the point cloud.
            transform: The coordinate transformation from the scene to the point cloud.
            subcollection: The name, or sequence of names, of the subcollection the
                point cloud is stored in. Defaults to ``'obsl'``.
            coordinate_space: Optional coordinate space for the point cloud. This will
                replace the existing coordinate space of the point cloud. Defaults to
                ``None``.

        Returns:
            The point cloud, opened for writing.

        Lifecycle: experimental
        """
        raise NotImplementedError()

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

    def get_transform_from_point_cloud(
        self, key: str, *, subcollection: str = "obsl"
    ) -> CoordinateTransform:
        """Returns the coordinate transformation from the requested point cloud to
        the scene.

        Args:
            key: The name of the point cloud.
            subcollection: The name, or sequence of names, of the subcollection the
                point cloud is stored in. Defaults to ``'obsl'``.

        Returns:
            Coordinate transform from the scene to the point cloud.

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
        raise NotImplementedError()

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
        raise NotImplementedError()

    def get_transform_to_point_cloud(
        self, key: str, *, subcollection: str = "obsl"
    ) -> CoordinateTransform:
        """Returns the coordinate transformation from the scene to a requested
        point cloud.

        Args:
            key: The name of the point cloud.
            subcollection: The name, or sequence of names, of the subcollection the
                point cloud is stored in. Defaults to ``'obsl'``.

        Returns:
            Coordinate transform from the scene to the requested point cloud.

        Lifecycle: experimental
        """
        raise NotImplementedError()
