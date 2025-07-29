# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""Column selections using index query condition or subarray. Internal use only."""

from __future__ import annotations

from typing import TYPE_CHECKING, Sequence

import numpy as np
import pyarrow as pa
from somacore import options

from . import pytiledbsoma as clib
from ._types import is_nonstringy_sequence
from ._util import _cast_domainish, pa_types_is_string_or_bytes, to_unix_ts

if TYPE_CHECKING:
    from ._soma_array import SOMAArray


class SOMAIndexValueFilter:
    _handle: clib.SOMAIndexColumnFilter
    _array: SOMAArray

    @classmethod
    def create(cls, array: SOMAArray, coords: options.SparseDFCoords | options.SparseNDCoords | options.DenseNDCoords):
        if not is_nonstringy_sequence(coords):
            raise TypeError(f"The coords type {type(coords)} must be a regular sequence, not str or bytes")

        if len(coords) > len(array._handle._handle.dimension_names):
            raise ValueError(
                f"The coords ({len(coords)} elements) must be shorter than ndim ({len(array._handle._handle.dimension_names)})",
            )

        qc = cls(array)
        for column_index, coord in enumerate(coords):
            qc.add_coordinate_selection(column_index, coord)

    def __init__(self, array: SOMAArray):
        self._array = array
        self._handle = clib.SOMAIndexQueryCondition(array)

    def add_coordinate_selection(
        self, column_index: int, coord: options.SparseDFCoord | options.SparseNDCoord | options.DenseNDCoord
    ) -> None:
        array_handle = self._array._handle._handle
        dim = array_handle.schema.field(column_index)
        domain = _cast_domainish(array_handle.domain())[column_index]
        if dim.metadata is not None and dim.metadata[b"dtype"].decode("utf-8") == "WKB":
            raise NotImplementedError("Support for adding a selection to a geometry column is not yet implemented.")
        if isinstance(coord, slice):
            if coord.step is not None and coord.step != 1:
                raise ValueError(
                    f"Invalid coordinate selection on column number '{column_index}'. Slice steps are not supported."
                )
            start = coord.start
            stop = coord.stop
            if pa_types_is_string_or_bytes(dim.type):
                dim_type = type(domain[0])
                # A ``None`` or empty start is always equivalent to empty str/bytes.
                start = coord.start or dim_type()
                if coord.stop is None:
                    # There's no way to specify "to infinity" for strings. We have to use this to get the nonempty domain
                    # and use that in the end.
                    ned = _cast_domainish(array_handle._handle.non_empty_domain())
                    _, stop = ned[column_index]
                else:
                    stop = coord.stop
                return self._handle.add_string_slice(start, stop)
            if pa.types.is_timestamp(dim.type) and (
                (start is None or isinstance(start, (np.datetime64, pa.TimestampScalar, int)))
                and (stop is None or isinstance(stop, (np.datetime64, pa.TimestampScalar, int)))
            ):
                # Convert datetime stamp to int64.
                istart = to_unix_ts(coord.start or domain[0])
                istop = to_unix_ts(coord.stop or domain[1])
                return self._handle.add_int64_slice(column_index, istart, istop)
            return self._handle.add_slice(column_index, start, stop)
        if isinstance(coord, (str, bytes)):
            return self._handle.add_string_points(column_index, coord)
        if isinstance(coord, (pa.Array, pa.ChunkedArray)):
            return self._handle.add_arrow_points(column_index, coord)
        if isinstance(coord, (Sequence, np.ndarray)):
            if isinstance(coord, np.ndarray) and coord.ndim != 1:
                raise ValueError(
                    f"Cannot set points on column index={column_index}. Can only use 1D numpy arrays, but got an array with {coord.ndim}."
                )
            self._handle.add_points(column_index, coord)
        raise TypeError("Unhandled type {dim.type} for index column named {dim.name}.")
