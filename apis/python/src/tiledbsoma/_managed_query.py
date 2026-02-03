# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""Managed query."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING, cast

import attrs
import numpy as np
import pyarrow as pa

import tiledbsoma.pytiledbsoma as clib

from . import _util
from ._core_options import PlatformConfig, SparseNDCoords
from ._types import is_nonstringy_sequence, is_slice_of

if TYPE_CHECKING:
    from ._soma_array import SOMAArray


@attrs.define(frozen=True)
class ManagedQuery:
    """Keep the lifetime of the SOMAArray tethered to ManagedQuery."""

    _array: SOMAArray
    _platform_config: PlatformConfig | None = None
    _handle: clib.ManagedQuery = attrs.field(init=False)
    _ref_store: list[object] = attrs.field(init=False, default=[])

    def __attrs_post_init__(self) -> None:
        array_handle = self._array._handle

        if self._platform_config is not None:
            cfg = array_handle.context().config()
            cfg.update(self._platform_config)
            ctx = clib.SOMAContext(cfg)
            object.__setattr__(self, "_handle", clib.ManagedQuery(array_handle, ctx))
        else:
            object.__setattr__(self, "_handle", clib.ManagedQuery(array_handle))

    def _set_coord_by_py_seq_or_np_array(self, dim: pa.Field, coord: object) -> None:
        if isinstance(coord, np.ndarray) and coord.ndim != 1:
            raise ValueError(f"only 1D numpy arrays may be used to index; got {coord.ndim}")

        column = self._array._handle.get_column(dim.name)

        try:
            set_dim_points = getattr(column, f"set_dim_points_{dim.type}")
        except AttributeError:
            # We have to handle this type specially below
            pass
        else:
            set_dim_points(self._handle, coord)
            return

        if _util.pa_types_is_string_or_bytes(dim.type):
            column.set_dim_points_string_or_bytes(self._handle, coord)
            return

        if pa.types.is_timestamp(dim.type):
            if not isinstance(coord, (tuple, list, np.ndarray)):
                raise ValueError(f"unhandled coord type {type(coord)} for index column named {dim.name}")

            icoord = [_util.to_unix_ts(e) for e in coord]
            column.set_dim_points_int64(self._handle, icoord)
            return

        raise ValueError(f"unhandled type {dim.type} for index column named {dim.name}")

    def _set_coord_by_numeric_slice(self, dim: pa.Field, dom: tuple[object, object], coord: slice) -> None:
        try:
            lo_hi = _util.slice_to_numeric_range(coord, dom)
        except _util.NonNumericDimensionError:
            return

        if not lo_hi:
            return

        column = self._array._handle.get_column(dim.name)

        try:
            set_dim_range = getattr(column, f"set_dim_ranges_{dim.type}")
            set_dim_range(self._handle, [lo_hi])
            return
        except AttributeError:
            return

    def set_coord(self, dim_idx: int, coord: object, axis_names: Sequence[str] | None = None) -> None:
        if coord is None:
            return

        array_handle = self._array._handle
        dim = array_handle.schema.field(dim_idx)
        dom = _util._cast_domainish(array_handle.domain())[dim_idx]

        column = array_handle.get_column(dim.name)

        if dim.metadata is not None and dim.metadata[b"dtype"].decode("utf-8") == "WKB":
            if axis_names is None:
                raise ValueError("Axis names are required to set geometry column coordinates")
            if not isinstance(dom[0], Mapping) or not isinstance(dom[1], Mapping):
                raise ValueError("Domain should be expressed per axis for geometry columns")

            self.set_geometry_coord(
                dim,
                cast("tuple[Mapping[str, float], Mapping[str, float]]", dom),
                coord,
                axis_names,
            )
            return

        if isinstance(coord, (str, bytes)):
            column.set_dim_points_string_or_bytes(self._handle, [coord])
            return

        if isinstance(coord, (pa.Array, pa.ChunkedArray)):
            column.set_dim_points_arrow(self._handle, coord)
            return

        if isinstance(coord, (Sequence, np.ndarray)):
            self._set_coord_by_py_seq_or_np_array(dim, coord)
            return

        if isinstance(coord, int):
            column.set_dim_points_int64(self._handle, [coord])
            return

        # Note: slice(None, None) matches the is_slice_of part, unless we also check
        # the dim-type part
        if (is_slice_of(coord, str) or is_slice_of(coord, bytes)) and _util.pa_types_is_string_or_bytes(dim.type):
            _util.validate_slice(coord)
            dim_type = type(dom[0])
            # A ``None`` or empty start is always equivalent to empty str/bytes.
            start = coord.start or dim_type()
            if coord.stop is None:
                # There's no way to specify "to infinity" for strings.
                # We have to get the nonempty domain and use that as the end.\
                ned = _util._cast_domainish(self._array._handle.non_empty_domain())
                _, stop = ned[dim_idx]
            else:
                stop = coord.stop
            column.set_dim_ranges_string_or_bytes(self._handle, [(start, stop)])
            return

        # Note: slice(None, None) matches the is_slice_of part, unless we also check
        # the dim-type part.
        if (
            is_slice_of(coord, np.datetime64) or is_slice_of(coord, pa.TimestampScalar) or is_slice_of(coord, int)
        ) and pa.types.is_timestamp(dim.type):
            _util.validate_slice(coord)

            # These timestamp types are stored in Arrow as well as TileDB as 64-bit
            # integers (with distinguishing metadata of course). For purposes of the
            # query logic they're just int64.
            istart = _util.to_unix_ts(coord.start or dom[0])
            istop = _util.to_unix_ts(coord.stop or dom[1])
            column.set_dim_ranges_int64(self._handle, [(istart, istop)])
            return

        if isinstance(coord, slice):
            _util.validate_slice(coord)
            if coord.start is None and coord.stop is None:
                return
            self._set_coord_by_numeric_slice(dim, dom, coord)
            return

        raise TypeError(f"unhandled type {dim.type} for index column named {dim.name}")

    def set_coords(self, coords: SparseNDCoords, axis_names: Sequence[str] | None = None) -> None:
        if not is_nonstringy_sequence(coords):
            raise TypeError(f"The coords type {type(coords)} must be a regular sequence, not str or bytes")

        if len(coords) > len(self._array._handle.dimension_names):
            raise ValueError(
                f"The coords ({len(coords)} elements) must be shorter than ndim"
                f" ({len(self._array._handle.dimension_names)})",
            )

        for i, coord in enumerate(coords):
            self.set_coord(i, coord, axis_names)

    def set_geometry_coord(
        self,
        dim: pa.Field,
        dom: tuple[Mapping[str, float], Mapping[str, float]],
        coord: object,
        axis_names: Sequence[str],
    ) -> None:
        if not isinstance(coord, Sequence):
            raise ValueError(f"unhandled coord type {type(coord)} for index column named {dim.name}")

        ordered_dom_min = [dom[0][axis] for axis in axis_names]
        ordered_dom_max = [dom[1][axis] for axis in axis_names]

        column = self._array._handle.get_column(dim.name)

        if all(is_slice_of(x, np.float64) for x in coord):
            range_min = []
            range_max = []
            for sub_coord, dom_min, dom_max in zip(coord, ordered_dom_min, ordered_dom_max):
                _util.validate_slice(sub_coord)
                lo_hi = _util.slice_to_numeric_range(sub_coord, (dom_min, dom_max))

                # None slices need to be set to the domain size
                if lo_hi is None:
                    range_min.append(dom_min)
                    range_max.append(dom_max)
                else:
                    range_min.append(lo_hi[0])
                    range_max.append(lo_hi[1])

            column.set_dim_ranges_double_array(self._handle, [(range_min, range_max)])
        elif all(isinstance(x, np.number) for x in coord):
            column.set_dim_points_double_array(self._handle, [coord])
        else:
            raise ValueError(f"Unsupported spatial coordinate type. Expected slice or float, found {type(coord)}")

    def set_column_data(self, dim_name: str, data: np.typing.NDArray) -> None:
        # store a reference to the data being written
        # libtiledbsoma will try not to copy any data to temporary buffers when writing data but the user is free
        # to pass temporary data objects to write. Submitting the write buffer to TileDB can happen after the
        # lifespan of the temporary object so ManagedQuery need to preserve a reference to each data object until
        # the query is submitted
        self._ref_store.append(data)

        self._handle.set_column_data(dim_name, data)

    def submit_batch(self, batch: pa.RecordBatch) -> None:
        # if array is remote we need to preserve the buffers until calling finalize
        if self._array._uri.startswith("tiledb://"):
            self._ref_store.append(batch)

        self._handle.submit_batch(batch)

    def submit_write(self) -> None:
        self._handle.submit_write()

        # if array is remote we need to preserve the buffers until calling finalize
        if not self._array._uri.startswith("tiledb://"):
            # clear stored data objects
            self._ref_store.clear()

    def submit_and_finalize(self) -> None:
        self._handle.submit_and_finalize()

        # clear stored data objects
        self._ref_store.clear()

    def finalize(self) -> None:
        self._handle.finalize()

        # clear stored data objects
        self._ref_store.clear()
