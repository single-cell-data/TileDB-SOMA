# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""Column selections using index query condition or subarray. Internal use only."""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING

import attrs
import numpy as np
import pyarrow as pa
from typing_extensions import Self

from . import pytiledbsoma as clib
from ._core_options import DenseCoord, DenseNDCoords, SparseDFCoord, SparseDFCoords, SparseNDCoord, SparseNDCoords
from ._types import is_nonstringy_sequence
from ._util import pa_types_is_string_or_bytes, to_unix_ts

if TYPE_CHECKING:
    from ._soma_array import SOMAArray


@attrs.define(frozen=True)
class CoordinateValueFilters:
    _array: SOMAArray
    _handle: clib.CoordinateValueFilters = attrs.field(init=False)

    @classmethod
    def create(cls, array: SOMAArray, coords: SparseDFCoords | SparseNDCoords | DenseNDCoords) -> Self:
        if not is_nonstringy_sequence(coords):
            raise TypeError(f"The coords type {type(coords)} must be a regular sequence, not str or bytes")

        if len(coords) > len(array._handle.dimension_names):
            raise ValueError(
                f"The coords ({len(coords)} elements) must be shorter than the number of index columns ({len(array._handle.dimension_names)})",
            )

        value_filter = cls(array)
        for column_index, coord in enumerate(coords):
            value_filter.add_coordinate_selection(column_index, coord)
        return value_filter

    def __attrs_post_init__(self) -> None:
        object.__setattr__(self, "_handle", clib.CoordinateValueFilters(self._array._handle))

    def add_coordinate_selection(self, column_index: int, coord: SparseDFCoord | SparseNDCoord | DenseCoord) -> None:
        array_handle = self._array._handle
        dim = array_handle.schema.field(column_index)
        if dim.metadata is not None and dim.metadata[b"dtype"].decode("utf-8") == "WKB":
            raise NotImplementedError("Support for adding a selection to a geometry column is not yet implemented.")

        if isinstance(coord, (pa.Array, pa.ChunkedArray)):
            self._handle.add_arrow_points(column_index, coord)
            return

        if isinstance(coord, slice):
            if coord.step is not None and coord.step != 1:
                raise ValueError(
                    f"Invalid coordinate selection on column number '{column_index}'. Slice steps are not supported."
                )

            if pa_types_is_string_or_bytes(dim.type):
                self._handle.add_slice_string(column_index, coord.start, coord.stop)
                return

            if pa.types.is_timestamp(dim.type):
                if (coord.start is None or isinstance(coord.start, (np.datetime64, pa.TimestampScalar, int))) and (
                    coord.stop is None or isinstance(coord.stop, (np.datetime64, pa.TimestampScalar, int))
                ):
                    # Convert datetime stamp to int64.
                    start = None if coord.start is None else to_unix_ts(coord.start)
                    stop = None if coord.stop is None else to_unix_ts(coord.stop)
                    self._handle.add_slice_int64(column_index, start, stop)
                    return
                raise TypeError("Unexpected type on column.")  # TODO Better error message before merging

            add_slice_function = getattr(self._handle, f"add_slice_{dim.type}")
            add_slice_function(column_index, coord.start, coord.stop)
            return

        if pa_types_is_string_or_bytes(dim.type):
            if isinstance(coord, np.ndarray) and coord.ndim != 1:
                raise ValueError(
                    f"Cannot set points on column index={column_index}. Can only use 1D numpy arrays, but got an array with {coord.ndim}."
                )
            self._handle.add_points_string(column_index, coord)
            return

        if isinstance(coord, (Sequence, np.ndarray)):
            if isinstance(coord, np.ndarray) and coord.ndim != 1:
                raise ValueError(
                    f"Cannot set points on column index={column_index}. Can only use 1D numpy arrays, but got an array with {coord.ndim}."
                )
            if pa.types.is_timestamp(dim.type):
                values = np.array(coord).astype(np.int64)
                self._handle.add_points_int64(column_index, values)
                return
            add_points_function = getattr(self._handle, f"add_points_{dim.type}")
            add_points_function(column_index, coord)
            return

        raise TypeError(f"Cannot add coordinate with type '{type(coord)}' to column {column_index}.")
