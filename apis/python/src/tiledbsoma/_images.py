# Copyright (c) 2024 TileDB, Inc,
#
# Licensed under the MIT License.
"""Implementation of a SOMA image collections."""

import abc
import json
from dataclasses import dataclass
from typing import Any, List, Optional, Sequence, Tuple, Union

import pyarrow as pa
import somacore
from somacore import coordinates, images, options
from typing_extensions import Final, Self

from . import _funcs, _tdb_handles
from ._collection import CollectionBase
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
        result_order: options.ResultOrderStr = somacore.ResultOrder.ROW_MAJOR,
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

    __slots__ = ("_axis_order", "_levels")
    _wrapper_type = _tdb_handles.Image2DCollectionWrapper

    _level_prefix: Final = "soma_level_"

    @dataclass
    class LevelProperties:
        """TODO: Add documentaiton for LevelProperties"""

        axis_order: Tuple[str, ...]
        name: str
        shape: Tuple[int, ...]

        def __init__(
            self, axis_order: Union[str, Sequence[str]], name: str, shape: Sequence[int]
        ):
            self.name = name
            self.axis_order = tuple(axis_order)
            self.shape = tuple(shape)
            self._channel_index = None
            for index, val in enumerate(self.axis_order):
                if val == "X":
                    self._x_index = index
                elif val == "Y":
                    self._y_index = index
                else:
                    # Must be channel.
                    self._channel_index = index

        def to_json(self) -> str:
            """Serialization of the class to JSON.

            Note that the level name is not serialzied with the other properties.
            """
            return json.dumps({"axis_order": self.axis_order, "shape": self.shape})

        @classmethod
        def from_json(cls, name: str, data: str) -> Self:
            """Construct class from a JSON string.

            Args:
                name: Name of the level.
                data: String storing JSON serialization of the level properties.
            """
            kwargs = json.loads(data)
            return cls(name=name, **kwargs)

        @property
        def height(self) -> int:
            """Number of pixels in the height of the image."""
            return self.shape[self._y_index]

        @property
        def width(self) -> int:
            """Number of pixels in the width of the image."""
            return self.shape[self._x_index]

    def __init__(
        self,
        handle: _tdb_handles.SOMAGroupWrapper[Any],
        **kwargs: Any,
    ):
        # Do generic SOMA collection initialization.
        super().__init__(handle, **kwargs)

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
        self._levels: List[Image2DCollection.LevelProperties] = []
        self._reset_levels()

    def _reset_levels(self) -> None:
        self._levels = [
            Image2DCollection.LevelProperties.from_json(
                key[len(self._level_prefix) :], val
            )
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
        props = Image2DCollection.LevelProperties(self.axis_order, key, shape)
        props_str = props.to_json()
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
        result_order: options.ResultOrderStr = somacore.ResultOrder.ROW_MAJOR,
        platform_config: Optional[options.PlatformConfig] = None,
    ) -> pa.Tensor:
        """TODO: Add read_image_level documentation"""
        raise NotImplementedError()
