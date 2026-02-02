# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.
"""Implementation of the SOMA image collection for spatial data."""

from __future__ import annotations

import abc
from collections.abc import MutableMapping, Sequence
from dataclasses import dataclass
from typing import Any, Final, Generic, TypeVar

import pyarrow as pa
from typing_extensions import Self

from . import base, coordinates, data, options

_DenseND = TypeVar("_DenseND", bound=data.DenseNDArray)
"""A particular implementation of a collection of DenseNDArrays."""
_RootSO = TypeVar("_RootSO", bound=base.SOMAObject)
"""The root SomaObject type of the implementation."""

_RO_AUTO = options.ResultOrder.AUTO
#
# Read types
#

_ReadData = TypeVar("_ReadData")
_UNBATCHED = options.BatchSize()


class PointCloudDataFrame(base.SOMAObject, metaclass=abc.ABCMeta):
    """A specialized SOMA DataFrame for storing collections of points in
    multi-dimensional space.

    The ``PointCloudDataFrame`` class is designed to efficiently store and query point
    data, where each point is represented by coordinates in one or more spatial
    dimensions (e.g., x, y, z) and may have additional columns for associated
    attributes.

    Lifecycle: experimental
    """

    __slots__ = ()
    soma_type: Final = "SOMAPointCloudDataFrame"  # type: ignore[misc]

    @classmethod
    @abc.abstractmethod
    def create(
        cls,
        uri: str,
        *,
        schema: pa.Schema,
        coordinate_space: Sequence[str] | coordinates.CoordinateSpace = (
            "x",
            "y",
        ),
        domain: Sequence[tuple[Any, Any] | None] | None = None,
        platform_config: options.PlatformConfig | None = None,
        context: Any | None = None,  # noqa: ANN401
    ) -> Self:
        """Creates a new ``PointCloudDataFrame`` at the given URI.

        The schema of the created point cloud  will include a column named
        ``soma_joinid`` of type ``pyarrow.int64``, with negative values disallowed, and
        at least one axis with numeric type.  If a ``soma_joinid`` column is present in
        the provided schema, it must be of the correct type.  If the ``soma_joinid``
        column is not provided, one will be added.


        The schema of the created point cloud must contain columns for the axes in the
        ``coordinate_space``. These columns will be index columns for the point cloud
        dataframe.

        Args:
            uri: The URI where the dataframe will be created.
            schema: Arrow schema defining the per-column schema. This schema
                must define all columns, including columns to be named as index
                columns.  If the schema includes types unsupported by the SOMA
                implementation, an error will be raised.
            coordinate_space: Either the coordinate space or the axis names for the
                coordinate space the point cloud is defined on.
            domain:
                An optional sequence of tuples specifying the domain of each
                index column. Each tuple must be a pair consisting of the
                minimum and maximum values storable in the index column.
                If provided, this sequence must have the same length as
                ``index_column_names``, and the index-column domain will be as
                specified.  If omitted entirely, or if ``None`` in a given
                dimension, the corresponding index-column domain will use an
                empty range, and data writes after that will fail with an
                exception.  Unless you have a particular reason not to, you
                should always provide the desired `domain` at create time: this
                is an optional but strongly recommended parameter.
            platform_config: platform-specific configuration; keys are SOMA
                implementation names.
            context: Other implementation-specific configuration.

        Returns:
            The newly created geometry dataframe, opened for writing.

        Lifecycle: experimental
        """
        raise NotImplementedError

    # Data operations

    @abc.abstractmethod
    def read(
        self,
        coords: options.SparseDFCoords = (),
        column_names: Sequence[str] | None = None,
        *,
        batch_size: options.BatchSize | None = _UNBATCHED,
        partitions: options.ReadPartitions | None = None,
        result_order: options.ResultOrderStr = _RO_AUTO,
        value_filter: str | None = None,
        platform_config: options.PlatformConfig | None = None,
    ) -> data.ReadIter[pa.Table]:
        """Reads a user-defined slice of data into Arrow tables.

        Args:
            coords: for each index dimension, which rows to read.
                Defaults to ``()``, meaning no constraint -- all IDs.
            column_names: the named columns to read and return.
                Defaults to ``None``, meaning no constraint -- all column names.
            batch_size: The size of batches that should be returned from a read.
                See :class:`options.BatchSize` for details.
            partitions: If present, specifies that this is part of
                a partitioned read, and which part of the data to include.
            result_order: the order to return results, specified as a
                :class:`~options.ResultOrder` or its string value.
            value_filter: an optional value filter to apply to the results.
                The default of ``None`` represents no filter. Value filter
                syntax is implementation-defined; see the documentation
                for the particular SOMA implementation for details.
            platform_config: platform-specific configuration; keys are SOMA
                implementation names.

        Returns:
            A :class:`ReadIter` of :class:`pa.Table`s.

        Lifecycle: experimental
        """
        raise NotImplementedError

    @abc.abstractmethod
    def read_spatial_region(
        self,
        region: options.SpatialRegion | None = None,
        column_names: Sequence[str] | None = None,
        *,
        region_transform: coordinates.CoordinateTransform | None = None,
        region_coord_space: coordinates.CoordinateSpace | None = None,
        batch_size: options.BatchSize = _UNBATCHED,
        partitions: options.ReadPartitions | None = None,
        result_order: options.ResultOrderStr = _RO_AUTO,
        value_filter: str | None = None,
        platform_config: options.PlatformConfig | None = None,
    ) -> SpatialRead[data.ReadIter[pa.Table]]:
        """Reads data intersecting an user-defined region of space into a
        :class:`SpatialRead` with data in Arrow tables.


        Args:
            region: The region to query. May be a box in the form
                [x_min, y_min, x_max, y_max] (for 2D images), a box in the form
                [x_min, y_min, z_min, x_max, y_max, z_max] (for 3D images), or
                a shapely Geometry.
            column_names: The named columns to read and return.
                Defaults to ``None``, meaning no constraint -- all column names.
            region_transform: An optional coordinate transform from the read region to the
                coordinate system of the spatial dataframe.
                Defaults to ``None``, meaning an identity transform.
            region_coord_space: An optional coordinate space for the region being read.
                Defaults to ``None``, coordinate space will be inferred from transform.
            batch_size: The size of batched reads.
                Defaults to `unbatched`.
            partitions: If present, specifies that this is part of a partitioned read,
                and which part of the data to include.
            result_order: the order to return results, specified as a
                :class:`~options.ResultOrder` or its string value.
            value_filter: an optional value filter to apply to the results.
                The default of ``None`` represents no filter. Value filter
                syntax is implementation-defined; see the documentation
                for the particular SOMA implementation for details.
            platform_config: platform-specific configuration; keys are SOMA
                implementation names.

        Returns:
            A :class:`SpatialRead` with :class:`ReadIter` of :class:`pa.Table`s data.

        Lifecycle: experimental
        """
        raise NotImplementedError

    @abc.abstractmethod
    def write(
        self,
        values: pa.RecordBatch | pa.Table,
        *,
        platform_config: options.PlatformConfig | None = None,
    ) -> Self:
        """Writes the data from an Arrow table to the persistent object.

        As duplicate index values are not allowed, index values already present
        in the object are overwritten and new index values are added.

        Args:
            values: An Arrow table containing all columns, including
                the index columns. The schema for the values must match
                the schema for the ``DataFrame``.
            platform_config: platform-specific configuration; keys are SOMA
                implementation names.

        Returns: ``self``, to enable method chaining.

        Lifecycle: experimental
        """
        raise NotImplementedError

    # Metadata operations

    @property
    @abc.abstractmethod
    def schema(self) -> pa.Schema:
        """The schema of the data in this dataframe.

        Lifecycle: experimental
        """
        raise NotImplementedError

    @property
    @abc.abstractmethod
    def index_column_names(self) -> tuple[str, ...]:
        """The names of the index (dimension) columns.

        Lifecycle: experimental
        """
        raise NotImplementedError

    @property
    @abc.abstractmethod
    def coordinate_space(self) -> coordinates.CoordinateSpace:
        """Coordinate space for this point cloud.

        Lifecycle: experimental
        """
        raise NotImplementedError

    @coordinate_space.setter
    @abc.abstractmethod
    def coordinate_space(self, value: coordinates.CoordinateSpace) -> None:
        """Coordinate space for this point cloud.

        Lifecycle: experimental
        """
        raise NotImplementedError

    @property
    @abc.abstractmethod
    def domain(self) -> tuple[tuple[Any, Any], ...]:
        """The allowable range of values in each index column.

        Returns: a tuple of minimum and maximum values, inclusive,
            storable on each index column of the dataframe.

        Lifecycle: experimental
        """
        raise NotImplementedError


class GeometryDataFrame(base.SOMAObject, metaclass=abc.ABCMeta):
    """A specialized SOMA object for storing complex geometries with spatial indexing.

    The ``GeometryDataFrame`` class is designed to store and manage geometric shapes such as
    polygons, lines, and multipoints, along with additional columns for associated attributes.

    Lifecycle: experimental
    """

    __slots__ = ()
    soma_type: Final = "SOMAGeometryDataFrame"  # type: ignore[misc]

    # Lifecycle

    @classmethod
    @abc.abstractmethod
    def create(
        cls,
        uri: str,
        *,
        schema: pa.Schema,
        coordinate_space: Sequence[str] | coordinates.CoordinateSpace = (
            "x",
            "y",
        ),
        domain: Sequence[tuple[Any, Any] | None] | None = None,
        platform_config: options.PlatformConfig | None = None,
        context: Any | None = None,  # noqa: ANN401
    ) -> Self:
        """Creates a new ``GeometryDataFrame`` at the given URI.

        The schema of the created geometry dataframe will include a column named
        ``soma_joinid`` of type ``pyarrow.int64``, with negative values
        disallowed, and a column named ``soma_geometry of type ``pyarrow.binary`` or
        ``pyarrow.large_binary``.  If a ``soma_joinid`` column or ``soma_geometry``
        are present in the provided schema, they must be of the correct type.  If
        either the ``soma_joinid`` column or ``soma_geometry`` column are not provided,
        one will be added.

        The geometry dataframe will be indexed using a spatial index for the
        ``soma_geometry`` column.

        Args:
            uri: The URI where the dataframe will be created.
            schema: Arrow schema defining the per-column schema. This schema
                must define all columns, including columns to be named as index
                columns.  If the schema includes types unsupported by the SOMA
                implementation, an error will be raised.
            coordinate_space: Either the coordinate space or the axis names for the
                coordinate space the point cloud is defined on.
            domain: An optional sequence of tuples specifying the domain of each
                index column. Two tuples must be provided for the ``soma_geometry``
                column which store the width followed by the height. Each tuple should
                be a pair consisting of the minimum and maximum values storable in the
                index column. If omitted entirely, or if ``None`` in a given dimension,
                the corresponding index-column domain will use the minimum and maximum
                possible values for the column's datatype.  This makes a dataframe
                growable.
            platform_config: platform-specific configuration; keys are SOMA
                implementation names.
            context: Other implementation-specific configuration.

        Returns:
            The newly created geometry dataframe, opened for writing.

        Lifecycle: experimental
        """
        raise NotImplementedError

    # Data operations

    @abc.abstractmethod
    def read(
        self,
        coords: options.SparseDFCoords = (),
        column_names: Sequence[str] | None = None,
        *,
        batch_size: options.BatchSize = _UNBATCHED,
        partitions: options.ReadPartitions | None = None,
        result_order: options.ResultOrderStr = _RO_AUTO,
        value_filter: str | None = None,
        platform_config: options.PlatformConfig | None = None,
    ) -> data.ReadIter[pa.Table]:
        """Reads a user-defined slice of data into Arrow tables.

        Args:
            coords: for each index dimension, which rows to read.
                Defaults to ``()``, meaning no constraint -- all IDs.
            column_names: the named columns to read and return.
                Defaults to ``None``, meaning no constraint -- all column names.
            batch_size: The size of batches that should be returned from a read.
                See :class:`options.BatchSize` for details.
            partitions: If present, specifies that this is part of
                a partitioned read, and which part of the data to include.
            result_order: the order to return results, specified as a
                :class:`~options.ResultOrder` or its string value.
            value_filter: an optional value filter to apply to the results.
                The default of ``None`` represents no filter. Value filter
                syntax is implementation-defined; see the documentation
                for the particular SOMA implementation for details.
            platform_config: platform-specific configuration; keys are SOMA
                implementation names.

        Returns:
            A :class:`ReadIter` of :class:`pa.Table`s.

        Lifecycle: experimental
        """
        raise NotImplementedError

    @abc.abstractmethod
    def read_spatial_region(
        self,
        region: options.SpatialRegion | None = None,
        column_names: Sequence[str] | None = None,
        *,
        region_transform: coordinates.CoordinateTransform | None = None,
        region_coord_space: coordinates.CoordinateSpace | None = None,
        batch_size: options.BatchSize = _UNBATCHED,
        partitions: options.ReadPartitions | None = None,
        result_order: options.ResultOrderStr = _RO_AUTO,
        value_filter: str | None = None,
        platform_config: options.PlatformConfig | None = None,
    ) -> SpatialRead[data.ReadIter[pa.Table]]:
        """Reads data intersecting an user-defined region of space into a
        :class:`SpatialRead` with data in Arrow tables.


        Args:
            region: The region to query. May be a box in the form
                [x_min, y_min, x_max, y_max] (for 2D images), a box in the form
                [x_min, y_min, z_min, x_max, y_max, z_max] (for 3D images), or
                a shapely Geometry.
            column_names: The named columns to read and return.
                Defaults to ``None``, meaning no constraint -- all column names.
            region_transform: An optional coordinate transform from the read region to the
                coordinate system of the spatial dataframe.
                Defaults to ``None``, meaning an identity transform.
            region_coord_space: An optional coordinate space for the region being read.
                Defaults to ``None``, coordinate space will be inferred from transform.
            batch_size: The size of batched reads.
                Defaults to `unbatched`.
            partitions: If present, specifies that this is part of a partitioned read,
                and which part of the data to include.
            result_order: the order to return results, specified as a
                :class:`~options.ResultOrder` or its string value.
            value_filter: an optional value filter to apply to the results.
                The default of ``None`` represents no filter. Value filter
                syntax is implementation-defined; see the documentation
                for the particular SOMA implementation for details.
            platform_config: platform-specific configuration; keys are SOMA
                implementation names.

        Returns:
            A :class:`SpatialRead` with :class:`ReadIter` of :class:`pa.Table`s data.

        Lifecycle: experimental
        """
        raise NotImplementedError

    @abc.abstractmethod
    def write(
        self,
        values: pa.RecordBatch | pa.Table,
        *,
        platform_config: options.PlatformConfig | None = None,
    ) -> Self:
        """Writes the data from an Arrow table to the persistent object.

        As duplicate index values are not allowed, index values already present
        in the object are overwritten and new index values are added.

        Args:
            values: An Arrow table containing all columns, including
                the index columns. The schema for the values must match
                the schema for the ``DataFrame``.
            platform_config: platform-specific configuration; keys are SOMA
                implementation names.

        Returns: ``self``, to enable method chaining.

        Lifecycle: experimental
        """
        raise NotImplementedError

    # Metadata operations

    @property
    @abc.abstractmethod
    def schema(self) -> pa.Schema:
        """The schema of the data in this dataframe.

        Lifecycle: experimental
        """
        raise NotImplementedError

    @property
    @abc.abstractmethod
    def index_column_names(self) -> tuple[str, ...]:
        """The names of the index (dimension) columns.

        Lifecycle: experimental
        """
        raise NotImplementedError

    @property
    @abc.abstractmethod
    def coordinate_space(self) -> coordinates.CoordinateSpace:
        """Coordinate space for this geometry dataframe.

        Lifecycle: experimental
        """
        raise NotImplementedError

    @coordinate_space.setter
    @abc.abstractmethod
    def coordinate_space(self, value: coordinates.CoordinateSpace) -> None:
        """Coordinate space for this geometry dataframe.

        Lifecycle: experimental
        """
        raise NotImplementedError

    @property
    @abc.abstractmethod
    def domain(self) -> tuple[tuple[Any, Any], ...]:
        """The allowable range of values in each index column.

        Returns: a tuple of minimum and maximum values, inclusive,
            storable on each index column of the dataframe.

        Lifecycle: experimental
        """
        raise NotImplementedError


class MultiscaleImage(
    base.SOMAObject,
    Generic[_DenseND, _RootSO],
    MutableMapping[str, _DenseND],
    metaclass=abc.ABCMeta,
):
    """A multiscale image with an extendable number of resolution levels.

    The multiscale image defines the top level properties. Each level must
    match the expected following properties:
    * number of channels
    * axis order
    * type

    Lifecycle: experimental
    """

    soma_type: Final = "SOMAMultiscaleImage"  # type: ignore[misc]
    __slots__ = ()

    # Lifecycle

    @classmethod
    @abc.abstractmethod
    def create(
        cls,
        uri: str,
        *,
        type: pa.DataType,
        level_shape: Sequence[int],
        level_key: str = "level0",
        level_uri: str | None = None,
        coordinate_space: Sequence[str] | coordinates.CoordinateSpace = (
            "x",
            "y",
        ),
        data_axis_order: Sequence[str] | None = None,
        platform_config: options.PlatformConfig | None = None,
        context: Any | None = None,  # noqa: ANN401
    ) -> Self:
        """Creates a new MultiscaleImage with one initial level.

        Args:
            uri: The URI where the collection will be created.
            type: The Arrow type to store the image data in the array.
                If the type is unsupported, an error will be raised.
            level_shape: The shape of the multiscale image for ``level=0``. Must
                include the channel dimension if there is one.
            level_key: The name for the ``level=0`` image. Defaults to ``level0``.
            level_uri: The URI for the ``level=0`` image. If the URI is an existing
                SOMADenseNDArray it must match have the shape provided by
                ``level_shape`` and type specified in ``type. If set to ``None``, the
                ``level_key`` will be used to construct a default child URI. For more
                on URIs see :meth:`collection.Collection.add_new_collection`.
            coordinate_space: Either the coordinate space or the axis names for the
                coordinate space the ``level=0`` image is defined on. This does not
                include the channel dimension, only spatial dimensions.
            data_axis_order: The order of the axes as stored on disk. Use
                ``soma_channel`` to specify the location of a channel axis. If no
                axis is provided, this defaults to the channel axis followed by the
                coordinate space axes in reverse order (e.g.
                ``("soma_channel", "y", "x")`` if ``coordinate_space=("x", "y")``).
            platform_config: platform-specific configuration; keys are SOMA
                implementation names.
            context: Other implementation-specific configuration.

        Returns:
            The newly created collection, opened for writing.

        Lifecycle: experimental
        """
        raise NotImplementedError

    @abc.abstractmethod
    def add_new_level(
        self,
        key: str,
        *,
        uri: str | None = None,
        shape: Sequence[int],
    ) -> _DenseND:
        """Add a new level in the multi-scale image.

        Parameters are as in :meth:`data.DenseNDArray.create`. The provided shape will
        be used to compute the scale between images and must correspond to the image
        size for the entire image. The image must be smaller than the ``level=0`` image.

        Lifecycle: experimental
        """
        raise NotImplementedError

    @abc.abstractmethod
    def set(
        self,
        key: str,
        value: _DenseND,
        *,
        use_relative_uri: bool | None = None,
    ) -> Self:
        """Sets a new level in the multi-scale image to be an existing SOMA
        :class:`data.DenseNDArray`.

        Args:
            key: The string key to set.
            value: The SOMA object to insert into the collection.
            use_relative_uri: Determines whether to store the collection
                entry with a relative URI (provided the storage engine
                supports it).
                If ``None`` (the default), will automatically determine whether
                to use an absolute or relative URI based on their relative
                location.
                If ``True``, will always use a relative URI. If the new child
                does not share a relative URI base, or use of relative URIs
                is not possible at all, the collection should raise an error.
                If ``False``, will always use an absolute URI.

        Returns: ``self``, to enable method chaining.

        Lifecycle: experimental
        """
        raise NotImplementedError

    # Data operations

    @abc.abstractmethod
    def read_spatial_region(
        self,
        level: int | str,
        region: options.SpatialRegion = (),
        *,
        channel_coords: options.DenseCoord = None,
        region_transform: coordinates.CoordinateTransform | None = None,
        region_coord_space: coordinates.CoordinateSpace | None = None,
        result_order: options.ResultOrderStr = _RO_AUTO,
        data_axis_order: Sequence[str] | None = None,
        platform_config: options.PlatformConfig | None = None,
    ) -> SpatialRead[pa.Tensor]:
        """Reads a user-defined region of space into a :class:`SpatialRead` with data
        in either an Arrow tensor or table.

        Reads the bounding box of the input region from the requested image level. This
        will return a :class:`SpatialRead` with the image data stored as a
        :class:`pa.Tensor`.

        Args:
            level: The image level to read the data from. May use index of the level
                or the image name.
            region: The region to query. May be a box in the form
                [x_min, y_min, x_max, y_max] (for 2D images), a box in the form
                [x_min, y_min, z_min, x_max, y_max, z_max] (for 3D images), or
                a shapely Geometry.
            channel_coords: An optional slice that defines the channel coordinates
                to read.
            region_transform: An optional coordinate transform that provides the
                transformation from the provided region to the reference level of this
                image. Defaults to ``None``.
            region_coord_space: An optional coordinate space for the region being read.
                The axis names must match the input axis names of the transform.
                Defaults to ``None``, coordinate space will be inferred from transform.
            data_axis_order: The order to return the data axes in. Use ``soma_channel``
                to specify the location of the channel coordinate.
            result_order: The order data to return results, specified as a
                :class:`~options.ResultOrder` or its string value. This is the result
                order the data is read from disk. It may be permuted if
                ``data_axis_order`` is not the default order.
            platform_config: platform-specific configuration; keys are SOMA
                implementation names.

        Returns:
            The data bounding the requested region as a :class:`SpatialRead` with
            :class:`pa.Tensor` data.
        """
        raise NotImplementedError

    # Metadata operations

    @property
    @abc.abstractmethod
    def coordinate_space(self) -> coordinates.CoordinateSpace:
        """Coordinate space for this multiscale image.

        Lifecycle: experimental
        """
        raise NotImplementedError

    @coordinate_space.setter
    @abc.abstractmethod
    def coordinate_space(self, value: coordinates.CoordinateSpace) -> None:
        """Coordinate space for this multiscale image.

        Lifecycle: experimental
        """
        raise NotImplementedError

    @property
    @abc.abstractmethod
    def data_axis_order(self) -> tuple[str, ...]:
        """The order of the axes for the images.

        Lifecycle: experimental
        """
        raise NotImplementedError

    @abc.abstractmethod
    def get_transform_from_level(self, level: int | str) -> coordinates.ScaleTransform:
        """Returns the transformation from user requested level to image reference level.

        Lifecycle: experimental
        """
        raise NotImplementedError

    @abc.abstractmethod
    def get_transform_to_level(self, level: int | str) -> coordinates.ScaleTransform:
        """Returns the transformation from the image reference level to the user
        requested level.

        Lifecycle: experimental
        """
        raise NotImplementedError

    @property
    @abc.abstractmethod
    def has_channel_axis(self) -> bool:
        """Returns if the images have an explicit channel axis.

        Lifecycle: experimental.
        """
        raise NotImplementedError

    @property
    @abc.abstractmethod
    def level_count(self) -> int:
        """The number of image levels stored in the MultiscaleImage.

        Lifecycle: experimental
        """
        raise NotImplementedError

    @abc.abstractmethod
    def level_shape(self, level: int | str) -> tuple[int, ...]:
        """The shape of the image at the specified level.

        Lifecycle: experimental
        """
        raise NotImplementedError

    @property
    @abc.abstractmethod
    def nchannels(self) -> int:
        """The number of channels.

        Lifecycle: experimental
        """
        raise NotImplementedError


@dataclass
class SpatialRead(Generic[_ReadData]):
    """Reader for spatial data.

    Args:
        data: The data accessor.
        data_coordinate_space: The coordinate space the read data is defined on.
        output_coordinate_space: The requested output coordinate space.
        coordinate_transform: A coordinate transform from the data coordinate space to
            the desired output coordinate space.

    Lifecycle: experimental
    """

    data: _ReadData
    data_coordinate_space: coordinates.CoordinateSpace
    output_coordinate_space: coordinates.CoordinateSpace
    coordinate_transform: coordinates.CoordinateTransform

    def __post_init__(self) -> None:
        if self.data_coordinate_space.axis_names != self.coordinate_transform.input_axes:
            raise ValueError("Input coordinate transform axis names do not match the data coordinate space.")
        if self.output_coordinate_space.axis_names != self.coordinate_transform.output_axes:
            raise ValueError("Output coordinate transform axis names do not match the output coordinate space.")
