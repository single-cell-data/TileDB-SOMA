#
# Copyright (c) 2024 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2024 TileDB, Inc.
#
# Licensed under the MIT License.
"""
Implementation of a SOMA MultiscaleImage.
"""

import json
import warnings
from typing import Any, Optional, Sequence, Tuple, Union

import attrs
import pyarrow as pa
import somacore
from somacore import (
    Axis,
    CoordinateSpace,
    CoordinateTransform,
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
class ImageProperties:
    """Properties for a single resolution level in a multiscale image.

    Lifecycle:
        Experimental.
    """

    name: str
    image_type: str
    shape: Tuple[int, ...] = attrs.field(converter=tuple)
    width: int = attrs.field(init=False)
    height: int = attrs.field(init=False)
    depth: Optional[int] = attrs.field(init=False)
    nchannels: Optional[int] = attrs.field(init=False)

    def __attrs_post_init__(self):  # type: ignore[no-untyped-def]
        if len(self.image_type) != len(set(self.image_type)):
            raise ValueError(
                f"Invalid image type '{self.image_type}'. Image type cannot contain "
                f"repeated values."
            )
        if len(self.image_type) != len(self.shape):
            raise ValueError(
                f"{len(self.image_type)} axis names must be provided for a multiscale "
                f"image with image type {self.image_type}."
            )

        nchannels: Optional[int] = None
        width: Optional[int] = None
        height: Optional[int] = None
        depth: Optional[int] = None
        for val, size in zip(self.image_type, self.shape):
            if val == "X":
                width = size
            elif val == "Y":
                height = size
            elif val == "Z":
                depth = size
            elif val == "C":
                nchannels = size
            else:
                raise SOMAError(f"Invalid image type '{self.image_type}'")
        if width is None or height is None:
            raise ValueError(
                f"Invalid image type '{self.image_type}'. Image type must include "
                f"'X' and 'Y'."
            )

        object.__setattr__(self, "nchannels", nchannels)
        object.__setattr__(self, "width", width)
        object.__setattr__(self, "height", height)
        object.__setattr__(self, "depth", depth)


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

    __slots__ = ("_schema", "_coord_space", "_levels")
    _wrapper_type = _tdb_handles.MultiscaleImageWrapper

    _level_prefix: Final = "soma_level_"

    # Lifecycle

    @classmethod
    def create(
        cls,
        uri: str,
        *,
        type: pa.DataType,
        reference_level_shape: Sequence[int],
        axis_names: Sequence[str] = ("c", "y", "x"),
        axis_types: Sequence[str] = ("channel", "height", "width"),
        platform_config: Optional[options.PlatformConfig] = None,
        context: Optional[SOMATileDBContext] = None,
        tiledb_timestamp: Optional[OpenTimestamp] = None,
    ) -> Self:
        """Creates a new ``MultiscaleImage`` at the given URI.

        Args:
            uri: The URI where the collection will be created.
            reference_level_shape: The shape of the reference level for the multiscale
                image. In most cases, this corresponds to the size of the image
                at ``level=0``.
            axis_names: The names of the axes of the image.
            axis_types: The types of the axes of the image. Must be the same length as
                ``axis_names``. Valid types are: ``channel``, ``height``, ``width``,
                and ``depth``.

        Returns:
            The newly created ``MultiscaleImage``, opened for writing.

        Lifecycle:
            Experimental.
        """
        # Warn about the experimental nature of the spatial classes.
        warnings.warn(SPATIAL_DISCLAIMER)

        context = _validate_soma_tiledb_context(context)
        if len(set(axis_types)) != len(axis_types):
            raise ValueError(
                "Invalid axis types {axis_types} - cannot have repeated values."
            )
        if len(axis_names) != len(axis_types):
            raise ValueError("Mismatched lengths for axis names and types.")
        axis_type_map = {"channel": "C", "height": "Y", "width": "X", "depth": "Z"}
        image_type = []
        for val in axis_types:
            try:
                image_type.append(axis_type_map[val])
            except KeyError as ke:
                raise ValueError(f"Invalid axis type name '{val}'.") from ke
        schema = MultiscaleImageSchema(
            ImageProperties(
                name="reference_level",
                image_type="".join(image_type),
                shape=tuple(reference_level_shape),  # type: ignore
            ),
            axis_names=tuple(axis_names),
            datatype=type,
        )

        # mypy false positive https://github.com/python/mypy/issues/5313
        coord_space = CoordinateSpace(
            tuple(Axis(name) for name in schema.get_coordinate_space_axis_names())  # type: ignore[misc]
        )
        schema_str = schema.to_json()
        coord_space_str = coordinate_space_to_json(coord_space)
        try:
            timestamp_ms = context._open_timestamp_ms(tiledb_timestamp)
            clib.SOMAGroup.create(
                uri=uri,
                soma_type=somacore.MultiscaleImage.soma_type,
                ctx=context.native_context,
                timestamp=(0, timestamp_ms),
            )
            handle = _tdb_handles.MultiscaleImageWrapper.open(
                uri, "w", context, tiledb_timestamp
            )
            handle.metadata[SOMA_MULTISCALE_IMAGE_SCHEMA] = schema_str
            handle.metadata[SOMA_COORDINATE_SPACE_METADATA_KEY] = coord_space_str
            return cls(
                handle,
                _dont_call_this_use_create_or_open_instead="tiledbsoma-internal-code",
            )
        except SOMAError as e:
            raise map_exception_for_create(e, uri) from None

    def __init__(
        self,
        handle: _tdb_handles.SOMAGroupWrapper[Any],
        **kwargs: Any,
    ):
        # Do generic SOMA collection initialization.
        super().__init__(handle, **kwargs)

        # Get schema for the multiscale image.
        try:
            schema_json = self.metadata[SOMA_MULTISCALE_IMAGE_SCHEMA]
        except KeyError as ke:
            raise SOMAError("Missing multiscale image schema metadata") from ke
        if isinstance(schema_json, bytes):
            schema_json = str(schema_json, "utf-8")
        if not isinstance(schema_json, str):
            raise SOMAError(
                f"Stored '{SOMA_MULTISCALE_IMAGE_SCHEMA}' metadata is unexpected "
                f"type {type(schema_json)!r}."
            )
        self._schema = MultiscaleImageSchema.from_json(schema_json)

        # Get the coordinate space.
        try:
            coord_space = self.metadata[SOMA_COORDINATE_SPACE_METADATA_KEY]
        except KeyError as ke:
            raise SOMAError("Missing coordinate space metadata") from ke
        self._coord_space = coordinate_space_from_json(coord_space)

        # Check schema and coordinate space have the same axis order
        schema_axes = self._schema.get_coordinate_space_axis_names()
        if schema_axes != self._coord_space.axis_names:
            raise SOMAError(
                f"Inconsistent axis names stored in metadata. Multiscale schema metadata"
                f" has coordinate axes '{schema_axes}', but the coordinate space "
                f"metadata has coordinate axes '{self._coord_space.axis_names}'"
            )

        # Get the image levels.
        # TODO: Optimize and push down to C++ level
        self._levels = [
            ImageProperties(name=key[len(self._level_prefix) :], **json.loads(val))
            for key, val in self.metadata.items()
            if key.startswith(self._level_prefix)
        ]
        self._levels.sort(key=lambda level: (-level.width, -level.height, level.name))

    @_funcs.forwards_kwargs_to(
        DenseNDArray.create, exclude=("context", "shape", "tiledb_timestamp")
    )
    def add_new_level(
        self,
        key: str,
        *,
        uri: Optional[str] = None,
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
        ref_props = self._schema.reference_level_properties
        shape = tuple(shape)
        if len(shape) != len(ref_props.shape):
            raise ValueError(
                f"New level must have {len(ref_props.shape)} dimensions, but shape "
                f"{shape} has {len(shape)} dimensions."
            )

        # Check, create, and store as metadata the new level image properties.
        props = ImageProperties(
            image_type=ref_props.image_type,
            name=key,
            shape=shape,  # type: ignore
        )
        if ref_props.nchannels is not None and ref_props.nchannels != props.nchannels:
            raise ValueError(
                f"New level must have {ref_props.nchannels}, but provided shape has "
                f"{props.nchannels} channels."
            )

        props_str = json.dumps(
            {
                "image_type": ref_props.image_type,
                "shape": shape,
            }
        )
        self.metadata[meta_key] = props_str

        # Add the level properties to level list.
        # Note: The names are guaranteed to be different from the earlier checks.
        for index, val in enumerate(self._levels):
            # Note: Name is unique, so guaranteed to be strict ordering.
            if (-props.width, -props.height, props.name) < (
                -val.width,
                -val.height,
                val.name,
            ):
                self._levels.insert(index, props)
                break
        else:
            self._levels.append(props)

        # Create and return new level array.

        return self._add_new_element(
            key,
            DenseNDArray,
            lambda create_uri: DenseNDArray.create(
                create_uri,
                context=self.context,
                tiledb_timestamp=self.tiledb_timestamp_ms,
                shape=props.shape,
                type=self._schema.datatype,
                **kwargs,
            ),
            uri,
        )

    # Data operations

    def read_spatial_region(
        self,
        level: Union[int, str],
        region: Optional[options.SpatialRegion] = None,
        *,
        channel_coords: options.DenseCoord = None,
        region_transform: Optional[CoordinateTransform] = None,
        region_coord_space: Optional[CoordinateSpace] = None,
        result_order: options.ResultOrderStr = options.ResultOrder.ROW_MAJOR,
        platform_config: Optional[options.PlatformConfig] = None,
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
            result_order: the order to return results, specified as a
                :class:`~options.ResultOrder` or its string value.

        Returns:
            The data bounding the requested region as a :class:`SpatialRead` with
            :class:`pa.Tensor` data.

        Lifecycle:
            Experimental.
        """
        # Get reference level. Check image is 2D.
        if self._schema.reference_level_properties.depth is not None:
            raise NotImplementedError(
                "Support for reading the levels of 3D images it not yet implemented."
            )

        # Check input query region type is supported.
        if (
            channel_coords is not None
            and self._schema.reference_level_properties.nchannels is None
        ):
            raise ValueError(
                "Invalide channel coordinate provided. This image has no channel "
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
                    f"Support for reading levels with a region tranform of type "
                    f"{type(region_transform)!r} is not yet supported."
                )
            # Create or check output coordinates.
            if region_coord_space is None:
                # mypy false positive https://github.com/python/mypy/issues/5313
                region_coord_space = CoordinateSpace(
                    tuple(Axis(axis_name) for axis_name in region_transform.input_axes)  # type: ignore[misc]
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
            self._schema.reference_level_properties.image_type,
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

    @property
    def axis_names(self) -> Tuple[str, ...]:
        """The name of the image axes.

        Lifecycle:
            Experimental.
        """
        return self._schema.axis_names

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

    def get_transform_from_level(self, level: Union[int, str]) -> ScaleTransform:
        """Returns the transformation from user requested level to the image reference
        level.

        Lifecycle:
            Experimental.
        """
        if isinstance(level, str):
            level_props = None
            for val in self._levels:
                if val.name == level:
                    level_props = val
                    break
            else:
                raise KeyError("No level with name '{level}'")
        else:
            level_props = self._levels[level]
        ref_level_props = self._schema.reference_level_properties
        if ref_level_props.depth is None:
            return ScaleTransform(
                input_axes=self._coord_space.axis_names,
                output_axes=self._coord_space.axis_names,
                scale_factors=[
                    ref_level_props.width / level_props.width,
                    ref_level_props.height / level_props.height,
                ],
            )
        assert level_props.depth is not None
        return ScaleTransform(
            input_axes=self._coord_space.axis_names,
            output_axes=self._coord_space.axis_names,
            scale_factors=[
                ref_level_props.width / level_props.width,
                ref_level_props.height / level_props.height,
                ref_level_props.depth / level_props.depth,
            ],
        )

    def get_transform_to_level(self, level: Union[int, str]) -> ScaleTransform:
        """Returns the transformation from the image reference level to the user
        requested level.

        Lifecycle:
            Experimental.
        """
        if isinstance(level, str):
            level_props = None
            for val in self._levels:
                if val.name == level:
                    level_props = val
                    break
            else:
                raise KeyError("No level with name '{level}'")
        else:
            level_props = self._levels[level]
        ref_level_props = self._schema.reference_level_properties
        if ref_level_props.depth is None:
            return ScaleTransform(
                input_axes=self._coord_space.axis_names,
                output_axes=self._coord_space.axis_names,
                scale_factors=[
                    level_props.width / ref_level_props.width,
                    level_props.height / ref_level_props.height,
                ],
            )
        assert level_props.depth is not None
        return ScaleTransform(
            input_axes=self._coord_space.axis_names,
            output_axes=self._coord_space.axis_names,
            scale_factors=[
                level_props.width / ref_level_props.width,
                level_props.height / ref_level_props.height,
                level_props.depth / ref_level_props.depth,
            ],
        )

    @property
    def image_type(self) -> str:
        """The order of the axes as stored in the data model.

        Lifecycle:
            Experimental.
        """
        return self._schema.reference_level_properties.image_type

    @property
    def level_count(self) -> int:
        """The number of image resolution levels stored in the ``MultiscaleImage``.

        Lifecycle:
            Experimental.
        """
        return len(self._levels)

    def level_properties(self, level: Union[int, str]) -> ImageProperties:
        """The properties of an image at the specified level.

        Lifecycle:
            Experimental.
        """
        if isinstance(level, str):
            raise NotImplementedError(
                "Support for getting level properties by name is not yet implemented."
            )  # TODO
        return self._levels[level]

    @property
    def reference_level(self) -> Optional[int]:
        """The index of image level that is used as a reference level.

        This will return ``None`` if no current image level matches the size of the
        reference level.

        Lifecycle:
            Experimental.
        """
        raise NotImplementedError()

    @property
    def reference_level_properties(self) -> "ImageProperties":
        """The image properties of the reference level.

        Lifecycle:
            Experimental.
        """
        return self._schema.reference_level_properties


# TODO: Push down to C++ layer
@attrs.define
class MultiscaleImageSchema:

    reference_level_properties: ImageProperties
    axis_names: Tuple[str, ...]
    datatype: pa.DataType

    def __attrs_post_init__(self):  # type: ignore[no-untyped-def]
        ndim = len(self.reference_level_properties.shape)
        if len(self.axis_names) != ndim:
            raise ValueError(
                f"Invalid axis names '{self.axis_names}'. {ndim} axis names must be "
                f"provided for a multiscale image with image type "
                f"{self.reference_level_properties.image_type}. "
            )

    # mypy false positive https://github.com/python/mypy/issues/5313
    def create_coordinate_space(self) -> CoordinateSpace:
        return CoordinateSpace(
            tuple(Axis(name) for name in self.get_coordinate_space_axis_names())  # type: ignore[misc]
        )

    def get_coordinate_space_axis_names(self) -> Tuple[str, ...]:
        # TODO: Setting axes and the coordinate space is going to be updated
        # in a future PR.
        x_name: Optional[str] = None
        y_name: Optional[str] = None
        z_name: Optional[str] = None
        for axis_name, axis_type in zip(
            self.axis_names, self.reference_level_properties.image_type
        ):
            if axis_type == "X":
                x_name = axis_name
            elif axis_type == "Y":
                y_name = axis_name
            elif axis_type == "Z":
                z_name = axis_name
        assert x_name is not None  # For mypy (already validated)
        assert y_name is not None  # For mypy (already validated)
        if z_name is None:
            return (x_name, y_name)
        else:
            return (x_name, y_name, z_name)

    def to_json(self) -> str:
        type_str = pyarrow_to_carrow_type(self.datatype)
        return json.dumps(
            {
                "name": self.reference_level_properties.name,
                "image_type": self.reference_level_properties.image_type,
                "shape": self.reference_level_properties.shape,
                "axis_names": self.axis_names,
                "datatype": type_str,
            }
        )

    @classmethod
    def from_json(cls, data: str) -> Self:
        kwargs = json.loads(data)
        axis_names = kwargs.pop("axis_names")
        type_str = kwargs.pop("datatype")
        type = carrow_type_to_pyarrow(type_str)
        return cls(ImageProperties(**kwargs), axis_names, type)
