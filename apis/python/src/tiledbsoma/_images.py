# Copyright (c) 2024 TileDB, Inc,
#
# Licensed under the MIT License.
"""Implementation of a SOMA image collections."""

import json
from dataclasses import dataclass
from typing import Any, Iterable, List, Optional, Tuple, Union

import pyarrow as pa
import somacore
from somacore import coordinates, images, options
from typing_extensions import Final, Self

from . import _funcs, _tdb_handles
from ._collection import CollectionBase
from ._dense_nd_array import DenseNDArray
from ._soma_object import AnySOMAObject
from ._types import OpenTimestamp
from .options import SOMATileDBContext


class Image2D(  # type: ignore[misc]  # __eq__ false positive
    CollectionBase[AnySOMAObject],
    images.Image2D[DenseNDArray, AnySOMAObject],
):
    """TODO: Add documentaiton for Image2D

    Lifecycle:
        Experimental.
    """

    __slots__ = "_levels"
    _wrapper_type = _tdb_handles.Image2DWrapper

    _level_prefix: Final = "soma_level_"

    @dataclass
    class LevelProperties:
        """TODO: Add documentaiton for LevelProperties"""

        axes: Tuple[str, ...]
        name: str
        shape: Tuple[int, ...]

        def __init__(
            self, axes: Union[str, Iterable[str]], name: str, shape: Iterable[int]
        ):
            self.axes = tuple(axes)
            self.shape = tuple(shape)
            if len(self.axes) != len(self.shape):
                raise ValueError()  # TODO: Add error message
            if set(self.axes) not in ({"X", "Y", "C"}, {"X", "Y"}):
                raise ValueError()  # TODO: Add error message
            self._channel_index = None
            for index, val in enumerate(self.axes):
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
            return json.dumps({"axes": self.axes, "shape": self.shape})

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
        super().__init__(handle, **kwargs)
        self._levels: List[Image2D.LevelProperties] = []

    def _reset_levels(self) -> None:
        self._levels = [
            Image2D.LevelProperties.from_json(key[len(self._level_prefix) :], val)
            for key, val in self.metadata.items()
            if key.startswith(self._level_prefix)
        ]
        self._levels.sort(key=lambda level: level.width, reverse=True)
        # TODO: Fix to sort by name if multiple values have the same width.

    @_funcs.forwards_kwargs_to(
        DenseNDArray.create, exclude=("context", "shape", "tiledb_timestamp")
    )
    def add_new_level(
        self,
        key: str,
        *,
        axes: Union[str, Iterable[str]],
        shape: Iterable[int],
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
        meta_key = f"soma_scale_{key}"
        if meta_key in self.metadata:
            raise KeyError(f"{key!r} already exists in {type(self)} scales")

        # Create the level property and store as metadata.
        props = Image2D.LevelProperties(axes, key, shape)
        props_str = props.to_json()
        self.metadata[meta_key] = props_str

        # Add the level properties to level list.
        # Note: The names are guaranteed to be different from the earlier checks.
        for index, val in enumerate(self._levels):
            if props.width > val.width or (
                props.width == val.width and props.name < val.name
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
    def level_count(self) -> int:
        return len(self._levels)

    def level_properties(self, level: int) -> images.Image2D.LevelProperties:
        return self._levels[level]

    @classmethod
    def open(
        cls,
        uri: str,
        mode: options.OpenMode = "r",
        *,
        tiledb_timestamp: Optional[OpenTimestamp] = None,
        context: Optional[SOMATileDBContext] = None,
        platform_config: Optional[options.PlatformConfig] = None,
        clib_type: Optional[str] = None,
    ) -> Self:
        """Opens this specific type of SOMA object."""
        obj = super().open(
            uri,
            mode,
            tiledb_timestamp=tiledb_timestamp,
            context=context,
            platform_config=platform_config,
            clib_type=clib_type,
        )
        obj._reset_levels()
        return obj

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
