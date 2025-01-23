# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.
"""
Implementation of a SOMA MultiscaleImage.
"""

from __future__ import annotations

import json
import warnings
from typing import Any, Dict, List, Sequence, Tuple, Union

import attrs
import pyarrow as pa
import somacore
from somacore import (
    CoordinateSpace,
    CoordinateTransform,
    IdentityTransform,
    ScaleTransform,
    options,
)
from typing_extensions import Final, Self

from . import _funcs, _tdb_handles
from . import pytiledbsoma as clib
from ._arrow_types import carrow_type_to_pyarrow, pyarrow_to_carrow_type
from ._constants import (
    SOMA_COORDINATE_SPACE_METADATA_KEY,
    SOMA_MULTISCALE_IMAGE_SCHEMA,
    SOMA_SPATIAL_VERSION_METADATA_KEY,
    SPATIAL_DISCLAIMER,
)
from ._dense_nd_array import DenseNDArray
from ._exception import SOMAError, map_exception_for_create
from ._soma_group import SOMAGroup
from ._soma_object import AnySOMAObject
from ._spatial_util import (
    coordinate_space_from_json,
    coordinate_space_to_json,
    process_image_region,
)
from ._types import OpenTimestamp
from .options import SOMATileDBContext
from .options._soma_tiledb_context import _validate_soma_tiledb_context


@attrs.define(frozen=True)
class _LevelProperties:
    """Properties for a single resolution level in a multiscale image."""

    name: str
    shape: Tuple[int, ...] = attrs.field(converter=tuple)


@attrs.define(frozen=True)
class _MultiscaleImageMetadata:
    """Helper class for reading/writing multiscale image metadata."""

    data_axis_permutation: Tuple[int, ...] = attrs.field(converter=tuple)
    has_channel_axis: bool
    shape: Tuple[int, ...] = attrs.field(converter=tuple)
    datatype: pa.DataType

    def to_json(self) -> str:
        type_str = pyarrow_to_carrow_type(self.datatype)
        return json.dumps(
            {
                "data_axis_permutation": self.data_axis_permutation,
                "shape": self.shape,
                "has_channel_axis": self.has_channel_axis,
                "datatype": type_str,
            }
        )

    @classmethod
    def from_json(cls, data: str) -> Self:
        kwargs = json.loads(data)
        type_str = kwargs.pop("datatype")
        type = carrow_type_to_pyarrow(type_str)
        return cls(datatype=type, **kwargs)


class MultiscaleImage(  # type: ignore[misc]  # __eq__ false positive
    SOMAGroup[DenseNDArray],
    somacore.MultiscaleImage[DenseNDArray, AnySOMAObject],
):
    """A multiscale image represented as a collection of images at multiple resolution levels.

    Each level of the multiscale image must have the following consistent properties:

    * **Number of Channels**: All levels must have the same number of channels.
    * **Axis Order**: The order of axes (e.g., channels, height, width) must be consistent across levels.

    Lifecycle:
        Experimental.
    """

    __slots__ = (
        "_coord_space",
        "_data_axis_permutation",
        "_datatype",
        "_has_channel_axis",
        "_levels",
    )
    _wrapper_type = _tdb_handles.MultiscaleImageWrapper

    _level_prefix: Final = "soma_level_"

    # Lifecycle

    @classmethod
    def create(
        cls,
        uri: str,
        *,
        type: pa.DataType,
        level_shape: Sequence[int],
        level_key: str = "level0",
        level_uri: str | None = None,
        coordinate_space: Union[Sequence[str], CoordinateSpace] = (
            "x",
            "y",
        ),
        data_axis_order: Sequence[str] | None = None,
        has_channel_axis: bool = True,
        platform_config: options.PlatformConfig | None = None,
        context: SOMATileDBContext | None = None,
        tiledb_timestamp: OpenTimestamp | None = None,
    ) -> Self:
        """Creates a new ``MultiscaleImage`` at the given URI.

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
            has_channel_axis: Save the image with a dedicated "channel" axis.
            platform_config: platform-specific configuration; keys are SOMA
                implementation names.
            context: Other implementation-specific configuration.
            tiledb_timestamp: set timestamp for created TileDB SOMA objects.


        Returns:
            The newly created ``MultiscaleImage``, with the initial image array open
            for writing.

        Lifecycle:
            Experimental.
        """
        # Warn about the experimental nature of the spatial classes.
        warnings.warn(SPATIAL_DISCLAIMER)

        context = _validate_soma_tiledb_context(context)

        # Create the coordinate space.
        if isinstance(coordinate_space, CoordinateSpace):
            axis_names = tuple(axis.name for axis in coordinate_space)
            axis_units = tuple(axis.unit for axis in coordinate_space)
        else:
            axis_names = tuple(coordinate_space)
            axis_units = tuple(len(axis_names) * [None])

        ndim = len(coordinate_space)
        if has_channel_axis:
            ndim += 1

        if len(level_shape) != ndim:
            channel_descript = (
                "with a channel axis" if has_channel_axis else "with no channel axis"
            )
            raise ValueError(
                f"Invalid shape {level_shape}. Expected {ndim} dimensions for a "
                f"multiscale image {channel_descript} on a coordinate space with "
                f"{len(coordinate_space)} dimensions."
            )
        if data_axis_order is None:
            axis_permutation = tuple(range(ndim - 1, -1, -1))
        else:
            axis_indices = {name: index for index, name in enumerate(axis_names)}
            if has_channel_axis:
                axis_indices["soma_channel"] = len(coordinate_space)
            if set(data_axis_order) != set(axis_indices.keys()):
                raise ValueError(
                    f"Invalid data axis order '{data_axis_order}'. Must be a "
                    f"permutation of the axes '{tuple(axis_indices.keys())}'."
                )
            axis_permutation = tuple(axis_indices[name] for name in data_axis_order)

        # The type ignore comments are to address a false positive in the attrs tuple
        # constructor.
        image_meta = _MultiscaleImageMetadata(
            data_axis_permutation=axis_permutation,  # type: ignore[arg-type]
            has_channel_axis=has_channel_axis,
            shape=level_shape,  # type: ignore[arg-type]
            datatype=type,
        )

        _image_meta_str = image_meta.to_json()
        try:
            timestamp_ms = context._open_timestamp_ms(tiledb_timestamp)
            clib.SOMAMultiscaleImage.create(
                uri=uri,
                axis_names=axis_names,
                axis_units=axis_units,
                ctx=context.native_context,
                timestamp=(0, timestamp_ms),
            )
            handle = _tdb_handles.MultiscaleImageWrapper.open(
                uri, "w", context, tiledb_timestamp
            )
            handle.metadata[SOMA_MULTISCALE_IMAGE_SCHEMA] = _image_meta_str
            multiscale = cls(
                handle,
                _dont_call_this_use_create_or_open_instead="tiledbsoma-internal-code",
            )
        except SOMAError as e:
            raise map_exception_for_create(e, uri) from None

        multiscale.add_new_level(
            level_key,
            uri=level_uri,
            shape=level_shape,
            platform_config=platform_config,
        )

        return multiscale

    def __init__(
        self,
        handle: _tdb_handles.SOMAGroupWrapper[Any],
        **kwargs: Any,
    ):
        # Do generic SOMA collection initialization.
        super().__init__(handle, **kwargs)

        try:
            spatial_encoding_version = self.metadata[SOMA_SPATIAL_VERSION_METADATA_KEY]
            if isinstance(spatial_encoding_version, bytes):
                spatial_encoding_version = str(spatial_encoding_version, "utf-8")
            if spatial_encoding_version not in {"0.1.0", "0.2.0"}:
                raise ValueError(
                    f"Unsupported MultiscaleImage with spatial encoding version "
                    f"{spatial_encoding_version}"
                )
        except KeyError as ke:
            raise SOMAError(
                "Missing spatial encoding version. May be deprecated experimental "
                "MultiscaleImage."
            ) from ke

        # Get the coordinate space.
        try:
            coord_space = self.metadata[SOMA_COORDINATE_SPACE_METADATA_KEY]
        except KeyError as ke:
            raise SOMAError("Missing coordinate space metadata") from ke
        self._coord_space = coordinate_space_from_json(coord_space)

        # Get the multiscale image specific metadata.
        try:
            metadata_json = self.metadata[SOMA_MULTISCALE_IMAGE_SCHEMA]
        except KeyError as ke:
            raise SOMAError("Missing multiscale image schema metadata") from ke
        if isinstance(metadata_json, bytes):
            metadata_json = str(metadata_json, "utf-8")
        if not isinstance(metadata_json, str):
            raise SOMAError(
                f"Stored '{SOMA_MULTISCALE_IMAGE_SCHEMA}' metadata is unexpected "
                f"type {type(metadata_json)!r}."
            )
        image_meta = _MultiscaleImageMetadata.from_json(metadata_json)
        self._data_axis_permutation = image_meta.data_axis_permutation
        self._has_channel_axis = image_meta.has_channel_axis
        self._datatype = image_meta.datatype

        # Get the image levels.
        # TODO: Optimize and push down to C++ level
        self._levels = [
            _LevelProperties(
                name=key[len(self._level_prefix) :], shape=tuple(json.loads(val))
            )
            for key, val in self.metadata.items()
            if key.startswith(self._level_prefix)
        ]
        self._levels.sort(
            key=lambda level: tuple(-val for val in level.shape) + (level.name,)
        )

    @_funcs.forwards_kwargs_to(
        DenseNDArray.create, exclude=("context", "shape", "tiledb_timestamp")
    )
    def add_new_level(
        self,
        key: str,
        *,
        uri: str | None = None,
        shape: Sequence[int],
        **kwargs: Any,
    ) -> DenseNDArray:
        """Adds a new resolution level to the ``MultiscaleImage``.

        Parameters are as in :meth:`DenseNDArray.create`. The provided shape will
        be used to compute the scale between images and must correspond to the image
        size for the entire image.

        Lifecycle:
            Experimental.
        """
        # Check if key already exists in either the collection or level metadata.
        if key in self:
            raise KeyError(f"{key!r} already exists in {type(self)}")
        meta_key = f"{self._level_prefix}{key}"
        if meta_key in self.metadata:
            raise KeyError(f"{key!r} already exists in {type(self)} scales")

        # Check if the shape is valid.
        shape = tuple(shape)
        ndim = len(self._data_axis_permutation)
        if len(shape) != ndim:
            raise ValueError(
                f"New level must have {ndim} dimensions, but shape {shape} has "
                f"{len(shape)} dimensions."
            )

        if self._has_channel_axis and len(self._levels) > 0:
            channel_index = self._data_axis_permutation.index(len(self._coord_space))
            expected_nchannels = self._levels[0].shape[channel_index]
            actual_nchannels = shape[channel_index]
            if actual_nchannels != expected_nchannels:
                raise ValueError(
                    f"New level must have {expected_nchannels}, but provided shape has "
                    f"{actual_nchannels} channels."
                )

        # Add the level properties to level list.
        # Note: The names are guaranteed to be different from the earlier checks.
        # Type ignore is for a false positive in the attrs tuple constructor.
        props = _LevelProperties(name=key, shape=shape)  # type: ignore[arg-type]
        for index, other in enumerate(self._levels):
            # Note: Name is unique, so guaranteed to be strict ordering.
            if tuple(-val for val in props.shape) + (props.name,) < tuple(
                -val for val in other.shape
            ) + (other.name,):
                self._levels.insert(index, props)
                break
        else:
            self._levels.append(props)

        shape_str = json.dumps(shape)
        self.metadata[meta_key] = shape_str

        # Create and return new level array.
        return self._add_new_element(
            key,
            DenseNDArray,
            lambda create_uri: DenseNDArray.create(
                create_uri,
                context=self.context,
                tiledb_timestamp=self.tiledb_timestamp_ms,
                shape=shape,
                type=self._datatype,
                **kwargs,
            ),
            uri,
        )

    def set(
        self,
        key: str,
        value: DenseNDArray,
        *,
        use_relative_uri: bool | None = None,
    ) -> Self:
        """Sets a new level in the multi-scale image to be an existing SOMA
        :class:`DenseNDArray`.

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
        raise NotImplementedError(
            "Support for setting external DenseNDArray objects to a MultiscaleImage "
            "is not yet implemented."
        )

    # Data operations

    def read_spatial_region(
        self,
        level: Union[int, str],
        region: options.SpatialRegion | None = None,
        *,
        channel_coords: options.DenseCoord = None,
        region_transform: CoordinateTransform | None = None,
        region_coord_space: CoordinateSpace | None = None,
        result_order: options.ResultOrderStr = options.ResultOrder.ROW_MAJOR,
        data_axis_order: Sequence[str] | None = None,
        platform_config: options.PlatformConfig | None = None,
    ) -> somacore.SpatialRead[pa.Tensor]:
        """Reads a user-defined spatial region from a specific level of the ``MultiscaleImage``.

        Retrieves the data within the specified region from the requested image level, returning
        a :class:`SpatialRead`, yielding :class:`pa.Tensor`s.

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

        Lifecycle:
            Experimental.
        """
        if data_axis_order is not None:
            raise NotImplementedError(
                "Support for altering the data axis order on read is not yet "
                "implemented."
            )

        # Get reference level. Check image is 2D.
        if len(self._coord_space) > 2:
            raise NotImplementedError(
                "Support for reading the levels of 3D images it not yet implemented."
            )

        # Check channel coords input is valid.
        if channel_coords is not None and not self._has_channel_axis:
            raise ValueError(
                "Invalid channel coordinate provided. This image has no channel "
                "dimension."
            )

        # Get the transformation for the group and the data coordinate space.
        # We may want to revisit copying the units for the data coordinate space.
        group_to_level = self.get_transform_to_level(level)
        data_coord_space = self.coordinate_space

        # Update transform and set region coordinate space.
        # - Add transformation from reference coord system to requested image level.
        # - Create or check the coordinate space for the input data region.
        if region_transform is None:
            if region_coord_space is not None:
                raise ValueError(
                    "Cannot specify the output coordinate space when region transform "
                    "is ``None``."
                )
            region_transform = group_to_level
            region_coord_space = data_coord_space
        else:
            if not isinstance(region_transform, ScaleTransform):
                raise NotImplementedError(
                    f"Support for reading levels with a region transform of type "
                    f"{type(region_transform)!r} is not yet supported."
                )
            # Create or check output coordinates.
            if region_coord_space is None:
                region_coord_space = CoordinateSpace.from_axis_names(
                    region_transform.input_axes
                )
            elif len(region_coord_space) != len(data_coord_space):
                raise ValueError(
                    "The number of output coordinates must match the number of "
                    "input coordinates."
                )
            if region_transform.output_axes != self._coord_space.axis_names:
                raise ValueError(
                    f"The output axes of '{region_transform.output_axes}' of the "
                    f"region transform must match the axes "
                    f"'{self._coord_space.axis_names}' of the coordinate space of "
                    f"this multiscale image."
                )
            region_transform = group_to_level @ region_transform
            assert isinstance(region_transform, ScaleTransform)

        # Convert coordinates to new coordinate system.
        coords, data_region, inv_transform = process_image_region(
            region,
            region_transform,
            channel_coords,
            self._data_axis_permutation,
        )

        # Get the array.
        array_name = level if isinstance(level, str) else self._levels[level].name
        try:
            array = self[array_name]
        except KeyError as ke:
            raise SOMAError(
                f"Unable to open the dense array with name '{array_name}'."
            ) from ke
        return somacore.SpatialRead(
            array.read(
                coords,
                result_order=result_order,
                platform_config=platform_config,
            ),
            data_coord_space,
            region_coord_space,
            inv_transform,
        )

    # Metadata operations
    def _level_properties(self, level: Union[int, str]) -> _LevelProperties:
        """The properties of an image at the specified level."""
        # by name
        # TODO could dynamically create a dictionary whenever a name-based
        # lookup is requested
        if isinstance(level, str):
            for val in self._levels:
                if val.name == level:
                    return val
            else:
                raise KeyError("No level with name '{level}'")

        # by index
        return self._levels[level]

    def _axis_order(self) -> List[int]:
        """Indices for accessing the data order for spatial axes."""
        axes = [
            index
            for index in range(len(self._data_axis_permutation))
            if self._data_axis_permutation[index] != len(self._coord_space)
        ]
        return sorted(axes, key=lambda index: self._data_axis_permutation[index])

    @property
    def coordinate_space(self) -> CoordinateSpace:
        """Coordinate space for this multiscale image.

        Lifecycle:
            Experimental.
        """
        return self._coord_space

    @coordinate_space.setter
    def coordinate_space(self, value: CoordinateSpace) -> None:
        """Coordinate space for this multiscale image.

        Lifecycle:
            Experimental.
        """
        if self._coord_space is not None:
            if value.axis_names != self._coord_space.axis_names:
                raise ValueError(
                    f"Cannot change axis names of a multiscale image. Existing axis "
                    f"names are {self._coord_space.axis_names}. New coordinate space "
                    f"has axis names {value.axis_names}."
                )
        self.metadata[SOMA_COORDINATE_SPACE_METADATA_KEY] = coordinate_space_to_json(
            value
        )
        self._coord_space = value

    @property
    def data_axis_order(self) -> Tuple[str, ...]:
        """The order of the axes for the resolution levels.

        Lifecycle:
            Experimental.
        """
        return tuple(
            (
                "soma_channel"
                if index == len(self._coord_space)
                else self._coord_space.axis_names[index]
            )
            for index in self._data_axis_permutation
        )

    def get_transform_from_level(self, level: Union[int, str]) -> ScaleTransform:
        """Returns the transformation from user requested level to the image reference
        level.

        Lifecycle:
            Experimental.
        """
        if level == 0 or level == self._level_properties(0).name:
            return IdentityTransform(
                input_axes=self._coord_space.axis_names,
                output_axes=self._coord_space.axis_names,
            )
        level_shape = self._level_properties(level).shape
        base_shape = self._levels[0].shape
        axis_indexer = self._axis_order()
        scale_factors = [
            base_shape[index] / level_shape[index] for index in axis_indexer
        ]
        return ScaleTransform(
            input_axes=self._coord_space.axis_names,
            output_axes=self._coord_space.axis_names,
            scale_factors=scale_factors,
        )

    def get_transform_to_level(self, level: Union[int, str]) -> ScaleTransform:
        """Returns the transformation from the image reference level to the user
        requested level.

        Lifecycle:
            Experimental.
        """
        if level == 0 or level == self._level_properties(0).name:
            return IdentityTransform(
                input_axes=self._coord_space.axis_names,
                output_axes=self._coord_space.axis_names,
            )
        level_shape = self._level_properties(level).shape
        base_shape = self._levels[0].shape
        axis_indexer = self._axis_order()
        scale_factors = [
            level_shape[index] / base_shape[index] for index in axis_indexer
        ]
        return ScaleTransform(
            input_axes=self._coord_space.axis_names,
            output_axes=self._coord_space.axis_names,
            scale_factors=scale_factors,
        )

    @property
    def has_channel_axis(self) -> bool:
        """Returns if the images have an explicit channel axis.

        Lifecycle:
            Experimental.
        """
        return self._has_channel_axis

    def levels(self) -> Dict[str, Tuple[str, Tuple[int, ...]]]:
        """Returns a mapping of {member_name: (uri, shape)}."""
        return {
            level.name: (self._contents[level.name].entry.uri, level.shape)
            for level in self._levels
        }

    @property
    def level_count(self) -> int:
        """The number of image resolution levels stored in the ``MultiscaleImage``.

        Lifecycle:
            Experimental.
        """
        return len(self._levels)

    def level_shape(self, level: Union[int, str]) -> Tuple[int, ...]:
        """The shape of the image at the specified resolution level.

        Lifecycle: experimental
        """
        if isinstance(level, str):
            for val in self._levels:
                if val.name == level:
                    return val.shape
            else:
                raise KeyError("No level with name '{level}'")

        # by index
        return self._levels[level].shape

    def level_uri(self, level: Union[int, str]) -> str:
        """The URI of the image at the specified resolution level.

        Lifecycle: experimental
        """
        if isinstance(level, int):
            level = self._levels[level].name
        return self._contents[level].entry.uri

    @property
    def nchannels(self) -> int:
        if self._has_channel_axis:
            index = self._data_axis_permutation.index(len(self._coord_space))
            return self._levels[0].shape[index]
        return 1
