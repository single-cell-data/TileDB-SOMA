# Copyright (c) 2024 TileDB, Inc,
#
# Licensed under the MIT License.
"""Implementation of a SOMA image collections."""

import abc
import json
from dataclasses import dataclass, field
from typing import Any, Optional, Sequence, Tuple, Union

import numpy as np
import pyarrow as pa
from somacore import ResultOrder, coordinates, images, options
from typing_extensions import Final

from . import _funcs, _tdb_handles
from ._collection import CollectionBase
from ._constants import SOMA_COORDINATE_SPACE_METADATA_KEY
from ._coordinates import (
    AffineCoordinateTransform,
    CoordinateSpace,
    CoordinateTransform,
)
from ._dense_nd_array import DenseNDArray
from ._exception import SOMAError
from ._soma_object import AnySOMAObject


class ImageCollection(  # type: ignore[misc]  # __eq__ false positive
    CollectionBase[AnySOMAObject],
    images.ImageCollection[DenseNDArray, AnySOMAObject],
    metaclass=abc.ABCMeta,
):

    @abc.abstractmethod
    def add_new_level(
        self,
        key: str,
        *,
        uri: Optional[str] = None,
        type: pa.DataType,
        shape: Sequence[int],
    ) -> DenseNDArray:
        """TODO: Add dcoumentation."""
        raise NotImplementedError()

    @property
    def axis_order(self) -> str:
        """The order of the axes in the stored images."""
        raise NotImplementedError()

    @property
    @abc.abstractmethod
    def level_count(self) -> int:
        """The number of image levels stored in the ImageCollection."""
        raise NotImplementedError()

    @abc.abstractmethod
    def level_properties(self, level: int) -> images.ImageCollection.LevelProperties:
        """The properties of an image at the specified level."""
        raise NotImplementedError()

    @abc.abstractmethod
    def read_level(
        self,
        level: int,
        coords: options.DenseNDCoords = (),
        *,
        transform: Optional[coordinates.CoordinateTransform] = None,
        result_order: options.ResultOrderStr = ResultOrder.ROW_MAJOR,
        platform_config: Optional[options.PlatformConfig] = None,
    ) -> pa.Tensor:
        """TODO: Add read_image_level documentation"""
        raise NotImplementedError()


class Image2DCollection(  # type: ignore[misc]  # __eq__ false positive
    ImageCollection, images.Image2DCollection
):
    """TODO: Add documentation for Image2DCollection

    Lifecycle:
        Experimental.
    """

    __slots__ = ("_axis_order", "_coord_space", "_levels", "_reference_shape")
    _wrapper_type = _tdb_handles.Image2DCollectionWrapper

    _level_prefix: Final = "soma_level_"

    @dataclass
    class LevelProperties:
        """TODO: Add documentaiton for LevelProperties"""

        name: str
        axis_order: str
        shape: Tuple[int, ...]
        width: int = field(init=False)
        height: int = field(init=False)
        nchannels: Optional[int] = field(init=False)

        def __post_init__(self):  # type: ignore[no-untyped-def]
            if len(self.axis_order) != len(self.shape):
                raise SOMAError()  # TODO Add error message
            self.nchannels = None
            for val, size in zip(self.axis_order, self.shape):
                if val == "X":
                    self.width = size
                elif val == "Y":
                    self.height = size
                elif val == "C":
                    self.nchannel = size
                else:
                    raise SOMAError(f"Invalid axis order '{self.axis_order}'")

    def __init__(
        self,
        handle: _tdb_handles.SOMAGroupWrapper[Any],
        **kwargs: Any,
    ):
        # Do generic SOMA collection initialization.
        super().__init__(handle, **kwargs)

        # Get the reference shape.
        ref_width = self.metadata.get("soma_image_reference_width")
        ref_height = self.metadata.get("soma_image_reference_height")
        if ref_width is not None and ref_height is not None:
            self._reference_shape: Optional[Tuple[int, ...]] = (
                int(ref_width),
                int(ref_height),
            )
            if any(val < 0 for val in self._reference_shape):
                raise SOMAError(
                    f"Decode invalid reference shape {self._reference_shape}"
                )
        elif not (ref_width is None and ref_height is None):
            raise SOMAError("Failed to fully decode reference shape.")
        else:
            self._reference_shape = None

        # Get the coordinate space.
        coord_space = self.metadata.get(SOMA_COORDINATE_SPACE_METADATA_KEY)
        if coord_space is None:
            self._coord_space: Optional[CoordinateSpace] = None
        else:
            self._coord_space = CoordinateSpace.from_json(coord_space)

        # Update the axis order.
        axis_order = self.metadata.get("soma_axis_order")
        if axis_order is None:
            self._axis_order: Optional[str] = None
        else:
            if isinstance(axis_order, bytes):
                axis_order = str(axis_order, "uft-8")
            if not isinstance(axis_order, str):
                raise SOMAError(
                    f"Stored Image2DCollection 'soma_axis_order' is unexpected type "
                    f"{type(axis_order)}."
                )
            if not (
                (len(axis_order) == 2 and set(axis_order) == {"X", "Y"})
                or (len(axis_order) == 3 and set(axis_order) == {"X", "Y", "C"})
            ):
                raise ValueError(
                    "Invalid axis order {axis_order}. The axis order must be a "
                    "permutation of 'CYX'."
                )
            self._axis_order = axis_order

        # Get the image levels.
        self._levels = [
            Image2DCollection.LevelProperties(name=key, **json.loads(val))
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
        shape: Sequence[int],
        uri: Optional[str] = None,
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

        # Create the level property and store as metadata.
        props = Image2DCollection.LevelProperties(
            axis_order=self.axis_order, name=key, shape=tuple(shape)
        )
        props_str = json.dumps({"axis_order": self.axis_order, "shape": shape})
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
                **kwargs,
            ),
            uri,
        )

    @property
    def axis_order(self) -> str:
        if self._axis_order is None:
            raise KeyError("The axis order is not set.")
        return self._axis_order

    @axis_order.setter
    def axis_order(self, value: str) -> None:
        if self._axis_order is not None and self._levels:
            raise ValueError(
                "The axis order is already set; if cannot be changed after adding "
                "images."
            )
        value = str(value).upper()
        if not (
            (len(value) == 2 and set(value) == {"X", "Y"})
            or (len(value) == 3 and set(value) == {"X", "Y", "C"})
        ):
            raise ValueError(
                "Invalid axis order {value}. The axis order must be a permutation of "
                "'CYX'."
            )
        self.metadata["soma_axis_order"] = value
        self._axis_order = value

    @property
    def coordinate_space(self) -> Optional[CoordinateSpace]:
        """Coordinate system for this image.

        The coordinate space is defined in order [width, height] even if the
        images are stored in a different order.
        """
        return self._coord_space

    @coordinate_space.setter
    def coordinate_space(self, value: CoordinateSpace) -> None:
        if self._reference_shape is None:
            raise SOMAError(
                "The reference shape must be set before the coordinate space"
            )
        if not isinstance(value, CoordinateSpace):
            raise TypeError(f"Invalid type {type(value).__name__}.")
        if len(value) != 2:
            raise ValueError("Coordinate space must have exactly 2 axes.")
        # TODO: Do we need some way to specify YX vs XY and propagate to
        # sub-images.
        self.metadata[SOMA_COORDINATE_SPACE_METADATA_KEY] = value.to_json()
        self._coord_space = value

    @property
    def level_count(self) -> int:
        return len(self._levels)

    def level_properties(self, level: int) -> images.ImageCollection.LevelProperties:
        return self._levels[level]

    def read_level(
        self,
        level: int,
        coords: options.DenseNDCoords = (),
        *,
        transform: Optional[coordinates.CoordinateTransform] = None,
        result_order: options.ResultOrderStr = ResultOrder.ROW_MAJOR,
        platform_config: Optional[options.PlatformConfig] = None,
    ) -> pa.Tensor:
        """TODO: Add read_image_level documentation"""
        raise NotImplementedError()

    @property
    def reference_shape(self) -> Tuple[int, ...]:
        """The shape of the reference level the coordinate system is defined on.

        The shape must be provide in order (width, height).
        """
        if self._reference_shape is None:
            raise SOMAError("Reference shape is not set.")
        return self._reference_shape

    @reference_shape.setter
    def reference_shape(self, value: Tuple[int, ...]) -> None:
        if self._reference_shape is not None:
            raise SOMAError("Cannot set reference shape; it is already set.")
        if len(value) != 2 or any(
            (not isinstance(val, int) or val < 0) for val in value
        ):
            raise ValueError("Invalid value")
        # TODO: Implement special handling of multi-value metadata
        self.metadata["soma_image_reference_width"] = value[0]
        self.metadata["soma_image_reference_height"] = value[1]
        self._reference_shape = (value[0], value[1])

    def transform_to_level(self, level: Union[str, int]) -> CoordinateTransform:
        if self._reference_shape is None or self._coord_space is None:
            raise SOMAError(
                "Cannot return transform if the reference shape and coordinate system "
                "for this image collection are not set."
            )
        if isinstance(level, str):
            for val in self._levels:
                if val.name == level:
                    level_props = val
                else:
                    raise KeyError("No level with name '{level}'")
        else:
            level_props = self._levels[level]
        width = level_props.width
        height = level_props.height
        # TODO: Add in a transformation just for scaling.
        # NOTE: Ignoring possibility of different axis order for different levels
        # since that will be removed.
        return AffineCoordinateTransform(
            input_axes=self._coord_space.axis_names,
            output_axes=self._coord_space.axis_names,
            matrix=np.array(
                [
                    [width / self._reference_shape[0], 0.0, 0.0],
                    [0.0, height / self._reference_shape[1], 0.0],
                    [0.0, 0.0, 1.0],
                ],
                dtype=np.float64,
            ),
        )
