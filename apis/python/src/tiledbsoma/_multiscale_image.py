# Copyright (c) 2024 TileDB, Inc,
#
# Licensed under the MIT License.
"""Implementation of a SOMA image collections."""

import json
from dataclasses import dataclass, field
from typing import Any, Optional, Sequence, Tuple, Union

import numpy as np
import pyarrow as pa
import shapely
import somacore
from somacore import (
    Axis,
    CoordinateSpace,
    CoordinateTransform,
    IdentityTransform,
    ResultOrder,
    ScaleTransform,
    options,
)
from typing_extensions import Final, Self

from . import _funcs, _tdb_handles
from . import pytiledbsoma as clib
from ._arrow_types import pyarrow_to_carrow_type
from ._constants import SOMA_COORDINATE_SPACE_METADATA_KEY, SOMA_MULTISCALE_IMAGE_SCHEMA
from ._coordinates import (
    coordinate_space_from_json,
    coordinate_space_to_json,
)
from ._dense_nd_array import DenseNDArray
from ._exception import SOMAError, map_exception_for_create
from ._soma_group import SOMAGroup
from ._soma_object import AnySOMAObject
from ._types import OpenTimestamp, is_nonstringy_sequence
from .options import SOMATileDBContext
from .options._soma_tiledb_context import _validate_soma_tiledb_context


@dataclass
class ImageProperties:
    """TODO: Add docstring"""

    name: str
    image_type: str
    shape: Tuple[int, ...]
    width: int = field(init=False)
    height: int = field(init=False)
    depth: Optional[int] = field(init=False)
    nchannels: Optional[int] = field(init=False)

    def __post_init__(self):  # type: ignore[no-untyped-def]
        if len(self.image_type) != len(set(self.image_type)):
            raise ValueError(
                f"Invalid image type '{self.image_type}'. Image type cannot contain "
                f"repeated values."
            )
        self.nchannels = None
        self.width = None  # type: ignore[assignment]
        self.height = None  # type: ignore[assignment]
        self.depth = None
        for val, size in zip(self.image_type, self.shape):
            if val == "X":
                self.width = size
            elif val == "Y":
                self.height = size
            elif val == "Z":
                self.depth = size
            elif val == "C":
                self.nchannels = size
            else:
                raise SOMAError(f"Invalid image type '{self.image_type}'")
        if self.width is None or self.height is None:
            raise ValueError(
                f"Invalid image type '{self.image_type}'. Image type must include "
                f"'X' and 'Y'."
            )
        if len(self.image_type) != len(self.shape):
            raise ValueError(
                f"{len(self.image_type)} axis names must be provided for a multiscale "
                f"image with image type {self.image_type}."
            )


class MultiscaleImage(  # type: ignore[misc]  # __eq__ false positive
    SOMAGroup[DenseNDArray],
    somacore.MultiscaleImage[DenseNDArray, AnySOMAObject],
):
    """TODO: Add documentation for MultiscaleImage

    Lifecycle:
        Experimental.
    """

    __slots__ = ("_schema", "_coord_space", "_levels")
    _wrapper_type = _tdb_handles.MultiscaleImageWrapper

    _level_prefix: Final = "soma_level_"

    @classmethod
    def create(
        cls,
        uri: str,
        *,
        type: pa.DataType,
        image_type: str = "CYX",
        reference_level_shape: Sequence[int],
        axis_names: Sequence[str] = ("c", "x", "y"),
        platform_config: Optional[options.PlatformConfig] = None,
        context: Optional[SOMATileDBContext] = None,
        tiledb_timestamp: Optional[OpenTimestamp] = None,
    ) -> Self:
        """TODO Add docstring

        Lifecycle:
            Experimental.
        """
        context = _validate_soma_tiledb_context(context)
        # TODO: Push down type to schema
        schema = MultiscaleImageSchema(
            ImageProperties(
                name="reference_level",
                image_type=image_type.upper(),
                shape=tuple(reference_level_shape),
            ),
            axis_names=tuple(axis_names),
            datatype=type,
        )

        coord_space = CoordinateSpace(
            tuple(Axis(name) for name in schema.get_coordinate_space_axis_names())
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
            ImageProperties(name=key, **json.loads(val))
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
        """Adds a new DenseNDArray to store the imagery for a new level


        TODO: explain how the parameters are used here. The remaining parameters
        are passed to the :meth:`DenseNDArray.create` method unchanged.
        """
        # Check if key already exists in either the collection or level metadata.
        if key in self:
            raise KeyError(f"{key!r} already exists in {type(self)}")
        meta_key = f"{self._level_prefix}{key}"
        if meta_key in self.metadata:
            raise KeyError(f"{key!r} already exists in {type(self)} scales")

        # Check if the shape is valid.
        ref_props = self._schema.reference_level_properties
        if len(shape) != len(tuple(shape)):
            raise ValueError(
                f"New level must have {len(shape)} dimensions, but shape {shape} has "
                f"{len(shape)} dimensions."
            )

        # Check, create, and store as metadata the new level image properties.
        props = ImageProperties(
            image_type=ref_props.image_type,
            name=key,
            shape=tuple(shape),
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

    @property
    def axis_names(self) -> Tuple[str, ...]:
        return self._schema.axis_names

    @property
    def coordinate_space(self) -> CoordinateSpace:
        """Coordinate system for this image.

        The coordinate space is defined in order [width, height] even if the
        images are stored in a different order.
        """
        return self._coord_space

    @coordinate_space.setter
    def coordinate_space(self, value: CoordinateSpace) -> None:
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

    def get_transformation_to_level(self, level: Union[str, int]) -> ScaleTransform:
        if isinstance(level, str):
            for val in self._levels:
                if val.name == level:
                    level_props = val
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

    def get_transformation_from_level(self, level: Union[str, int]) -> ScaleTransform:
        if isinstance(level, str):
            for val in self._levels:
                if val.name == level:
                    level_props = val
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

    @property
    def image_type(self) -> str:
        return self._schema.reference_level_properties.image_type

    @property
    def level_count(self) -> int:
        return len(self._levels)

    def level_properties(self, level: Union[int, str]) -> somacore.ImageProperties:
        if isinstance(level, str):
            raise NotImplementedError(
                "Support for getting level properties by name is not yet implemented."
            )  # TODO
        return self._levels[level]

    def read_level(
        self,
        level: Union[int, str],
        region: options.ImageCoords = (),
        *,
        transform: Optional[CoordinateTransform] = None,
        region_coord_space: Optional[CoordinateSpace] = None,
        apply_mask: bool = False,
        result_order: options.ResultOrderStr = ResultOrder.ROW_MAJOR,
        platform_config: Optional[options.PlatformConfig] = None,
    ) -> somacore.SpatialRead[pa.Tensor]:
        """Reads a user-defined dense slice of the array.

        If provided, the transform is interpreted as the coordinate transformation
        from the reference multiscale level to the requested region.
        """
        # Applying a mask in not yet supported.
        if apply_mask:
            raise NotImplementedError(
                "Support for applying a mask to the image is not yet implemented."
            )

        # Check input query region type is supported.
        if isinstance(region, shapely.GeometryType):
            raise NotImplementedError(
                "Support for querying by geometry is not yet implemented."
            )
        if not is_nonstringy_sequence(region):
            raise TypeError(
                f"non-geometry region type {type(region)} must be a regular sequence,"
                " not str or bytes"
            )

        # Get reference level. Check image is 2D.
        ref_level_props = self._schema.reference_level_properties
        if ref_level_props.depth is not None:
            raise NotImplementedError(
                "Support for reading the levels of 3D images it not yet implemented."
            )

        # Get the level properties for scaling.
        if isinstance(level, str):
            raise NotImplementedError("Accessing level by name is not yet supported.")
        level_props = self._levels[level]

        # Get the data coordinate space
        def scaled_axis(
            axis: Axis, a: Union[int, np.float64], b: Union[int, np.float64]
        ) -> Axis:
            return Axis(
                axis.name,
                axis.units,
                None if axis.scale is None else a * axis.scale / b,
            )

        x_axis = self.coordinate_space.axes[0]
        y_axis = self.coordinate_space.axes[1]
        data_coord_space = CoordinateSpace(
            (
                scaled_axis(x_axis, level_props.width, ref_level_props.width),
                scaled_axis(y_axis, level_props.height, ref_level_props.height),
            )
        )

        # Get the transform from data coord to multiscale image reference level
        if transform is None:
            if region_coord_space is not None:
                raise ValueError(
                    "Cannot specify the output coordinate space when transform is "
                    "``None``."
                )
            transform = IdentityTransform(
                data_coord_space.axis_names, data_coord_space.axis_names
            )
            region_coord_space = data_coord_space
        else:
            if not isinstance(transform, ScaleTransform):
                raise NotImplementedError(
                    f"Support for reading levels with a tranform of type "
                    f"{type(transform)!r} is not yet supported."
                )
            # Create or check output coordinates.
            if region_coord_space is None:
                region_coord_space = data_coord_space  # TODO: Fix this
            elif len(region_coord_space) != len(data_coord_space):
                raise ValueError(
                    "The number of output coordinates must match the number of "
                    "input coordinates."
                )

        # TODO: Quick and dirty implementation - clean this up later.
        def transform_coord(
            axis_type: str, coord: options.DenseCoord
        ) -> options.DenseCoord:
            if coord is None or axis_type == "C":
                return coord
            if axis_type == "X":
                axis_index = 0
                scale = level_props.width / ref_level_props.width
            elif axis_type == "Y":
                axis_index = 1
                scale = level_props.height / ref_level_props.width
            else:
                raise ValueError(f"unexpected axis type '{axis_type}'")
            scale *= (
                transform.scale_factors
                if transform.isotropic
                else transform.scale_factors[axis_index]  # type: ignore[index]
            )

            if isinstance(coord, int):
                return int(np.round(scale * coord))
            if isinstance(coord, slice):
                if coord.step is not None:
                    raise NotImplementedError(
                        "Support for steps in image read is not supported."
                    )
                return slice(
                    None if coord.start is None else int(np.floor(scale * coord.start)),
                    None if coord.stop is None else int(np.ceil(scale * coord.stop)),
                )
            raise TypeError(f"Unexpected type coordinate {type(coord)!r} in region")

        # transform = (data_to_reference + offset) * transform
        # TODO: Apply inverse transform on requested region.
        # Maybe use shapely?
        # Note: hard-coded assuming transform only supports ScaleTransfrom
        if not isinstance(transform, IdentityTransform):
            coords: options.DenseNDCoords = tuple(
                transform_coord(axis_type, coord)
                for axis_type, coord in zip(level_props.image_type, region)
            )

        # Get the array.
        try:
            array = self[level_props.name]
        except KeyError as ke:
            raise SOMAError(f"Unable to open the dense array at level {level}.") from ke
        return somacore.SpatialRead(
            array.read(
                coords,
                result_order=result_order,
                platform_config=platform_config,
            ),
            data_coord_space,
            region_coord_space,
            transform,
        )

    @property
    def reference_level(self) -> Optional[int]:
        """Returns the level of the image the coordinate system is defined with
        respect to."""
        raise NotImplementedError()

    @property
    def reference_level_properties(self) -> somacore.ImageProperties:
        """The shape of the reference level the coordinate system is defined on.

        The shape must be provide in order (width, height).
        """
        return self._schema.reference_level_properties


# TODO: Push down to C++ layer
@dataclass
class MultiscaleImageSchema:

    reference_level_properties: ImageProperties
    axis_names: Tuple[str, ...]
    datatype: pa.DataType

    def __post_init__(self):  # type: ignore[no-untyped-def]
        ndim = len(self.reference_level_properties.shape)
        if len(self.axis_names) != ndim:
            raise ValueError(
                f"Invalid axis names '{self.axis_names}'. {ndim} axis names must be "
                f"provided for a multiscale image with image type "
                f"{self.reference_level_properties.image_type}. "
            )

    def create_coordinate_space(self) -> CoordinateSpace:
        return CoordinateSpace(
            tuple(Axis(name) for name in self.get_coordinate_space_axis_names())
        )

    def get_coordinate_space_axis_names(self) -> Tuple[str, ...]:
        channel_index = self.reference_level_properties.image_type.find("C")
        if channel_index == -1:
            return self.axis_names
        return tuple(
            self.axis_names[:channel_index] + self.axis_names[channel_index + 1 :]
        )

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
        _carrow_to_pyarrow = {
            "c": pa.int8(),
            "s": pa.int16(),
            "i": pa.int32(),
            "l": pa.int64(),
            "C": pa.uint8(),
            "S": pa.uint16(),
            "I": pa.uint32(),
            "L": pa.uint64(),
            "f": pa.float32(),
            "g": pa.float64(),
        }
        type = _carrow_to_pyarrow[type_str]
        return cls(ImageProperties(**kwargs), axis_names, type)
