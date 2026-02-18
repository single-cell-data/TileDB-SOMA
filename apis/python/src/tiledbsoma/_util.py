# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

from __future__ import annotations

import datetime
import pathlib
import time
import urllib.parse
from concurrent.futures import Future
from itertools import zip_longest
from string import ascii_lowercase, ascii_uppercase, digits
from typing import Any, TypeVar

import numpy as np
import pandas as pd
import pyarrow as pa

from . import pytiledbsoma as clib
from ._core_options import DenseCoord, DenseNDCoords, IOfN, ReadPartitions, ResultOrder, ResultOrderStr
from ._types import DataProtocol, OpenTimestamp, Slice, is_slice_of


def get_start_stamp() -> float:
    """Returns information about start time of an event.

    Nominally float seconds since the epoch, but articulated here
    as being compatible with the format_elapsed function.
    """
    return time.time()


def format_elapsed(start_stamp: float, message: str) -> str:
    """Returns the message along with an elapsed-time indicator,
    with end time relative to start start from ``get_start_stamp``.

    Used for annotating elapsed time of a task.
    """
    return f"{message} TIME {time.time() - start_stamp:.3f} seconds"


def is_local_path(path: str) -> bool:
    if path.startswith("file://"):
        return True
    return "://" not in path


def make_relative_path(uri: str, relative_to: str) -> str:
    """Returns a URI relative to another URI. If not possible, raise a ValueError.

    This function assumes that the URI scheme follows posix path conventions
    and only contains a scheme, netloc and path. It does not handle query params,
    fragments, etc.

    The default scheme, if one is not specified, is assumed to be `file`.
    """
    p_uri = urllib.parse.urlparse(uri)
    p_relative_to = urllib.parse.urlparse(relative_to)

    uri_scheme = p_uri.scheme if p_uri.scheme else "file"
    relative_to_scheme = p_relative_to.scheme if p_relative_to.scheme else "file"
    if uri_scheme != relative_to_scheme:
        raise ValueError("Unable to make relative path between URIs with different scheme")

    return pathlib.PurePath(p_uri.path).relative_to(p_relative_to.path).as_posix()


def is_relative_uri(uri: str) -> bool:
    """Detects whether a provided URI is a child relative URI to the parent."""
    return "://" not in uri and not uri.startswith("/")


def uri_joinpath(base: str, path: str) -> str:
    """Join a path to a URI.

    Supports relative paths for ``file`` or unspecified schemes, assuming
    they are file system paths.  Assumes NO support for relative paths
    otherwise.
    """
    p_base = urllib.parse.urlparse(base)
    parts = [*p_base]

    if len(path) == 0:
        return base

    if not p_base.scheme or p_base.scheme == "file":
        # if a file path, just use pathlib. This is significantly more
        # permissive than it should be, given that `file://` URIs are
        # only absolute.
        parts[2] = pathlib.PurePath(p_base.path).joinpath(path).as_posix()
    else:
        if ".." in path:
            raise ValueError("Relative paths unsupported")
        if path.startswith("/"):
            # if absolute, just use the path
            parts[2] = path
        else:
            # join, being careful about extraneous path sep
            if parts[2].endswith("/"):
                parts[2] = parts[2] + path
            else:
                parts[2] = parts[2] + "/" + path

    return urllib.parse.urlunparse(parts)


def validate_slice(slc: Slice[Any]) -> None:
    """Checks that a slice has no step and is not inverted."""
    if slc.step is not None:
        raise ValueError("slice steps are not supported")
    if slc.start is None or slc.stop is None:
        # All half-specified slices are valid.
        return

    if isinstance(slc.stop, pa.TimestampScalar) or isinstance(slc.start, pa.TimestampScalar):
        if to_unix_ts(slc.stop) < to_unix_ts(slc.start):
            raise ValueError(f"slice start ({slc.start!r}) must be <= slice stop ({slc.stop!r})")
        return

    if slc.stop < slc.start:
        raise ValueError(f"slice start ({slc.start!r}) must be <= slice stop ({slc.stop!r})")


_T = TypeVar("_T")


class NonNumericDimensionError(TypeError):
    """Raised when trying to get a numeric range for a non-numeric dimension."""


def slice_to_numeric_range(slc: Slice[Any], domain: tuple[_T, _T]) -> tuple[_T, _T] | None:
    """Constrains the given slice to the ``domain`` for numeric dimensions.

    We assume the slice has already been validated by validate_slice.
    """
    if slc == slice(None):
        return None

    domain_start, domain_stop = domain
    if isinstance(domain_stop, (str, bytes)):
        # Strings don't have a real "domain" so we can't handle them
        # the same way that we handle numeric types.
        raise NonNumericDimensionError("only numeric dimensions supported")

    # TODO: with future C++ improvements, move half-slice logic to SOMAArrayReader
    start = domain_start if slc.start is None else max(slc.start, domain_start)
    stop = domain_stop if slc.stop is None else min(slc.stop, domain_stop)

    # Lint says the left-hand and right-hand sides are both unions
    if stop < start:  # type: ignore[operator]
        # With the above, we have guaranteed that at least one bound will
        # include the domain.  If we get here, that means that the other bound
        # never included it (e.g. stop == slc.stop < domain_start == start).
        raise ValueError(f"slice [{slc.start!r}:{slc.stop!r}] does not overlap [{domain_start!r}:{domain_stop!r}]")

    return start, stop


def dense_indices_to_shape(
    coords: DenseNDCoords,
    array_shape: tuple[int, ...],
    result_order: ResultOrder,
) -> tuple[int, ...]:
    """Given a subarray index specified as a tuple of per-dimension slices or scalars
    (e.g., ``([:], 1, [1:2])``), and the shape of the array, return the shape of
    the subarray. Note that the number of coordinates may be less than or equal
    to the number of dimensions in the array.
    """
    if len(coords) > len(array_shape):
        raise ValueError(f"coordinate length ({len(coords)}) must be <= array dimension count ({len(array_shape)})")

    shape = tuple(dense_index_to_shape(coord, extent) for coord, extent in zip_longest(coords, array_shape))
    if result_order == ResultOrder.ROW_MAJOR:
        return shape
    return tuple(reversed(shape))


def dense_index_to_shape(coord: DenseCoord, array_length: int) -> int:
    """Given a subarray per-dimension index specified as a slice or scalar (e.g, ``[:], 1, [1:2]``),
    and the shape of the array in that dimension, return the shape of the subarray in
    that dimension.

    Note that Python slice semantics are right-endpoint-exclusive whereas SOMA slice semantics are
    doubly inclusive.
    """
    if coord is None:
        return array_length
    if isinstance(coord, int):
        return 1
    if is_slice_of(coord, int):
        # We verify that ``step`` is None elsewhere, so we can always assume
        # that we're asked for a continuous slice.
        stop = array_length if coord.stop is None else min(coord.stop + 1, array_length)
        return stop - (coord.start or 0)

    raise TypeError(f"coordinate {coord} must be integer or integer slice")


def check_type(
    name: str,
    actual_value: Any,  # noqa: ANN401
    expected_types: tuple[type[Any], ...],
) -> None:
    """Verifies the type of an argument, or produces a useful error message."""
    if not isinstance(actual_value, expected_types):
        if len(expected_types) == 1:
            raise TypeError(f"expected {name} argument to be of type {expected_types[0]}; got {type(actual_value)}")
        raise TypeError(f"expected {name} argument to be one of {expected_types!r}; got {type(actual_value)}")


def check_unpartitioned(partitions: ReadPartitions | None) -> None:
    """Ensures that we're not being asked for a partitioned read.

    Because we currently don't support partitioned reads, we should reject all
    reads that request partitions to avoid giving the user duplicate data across
    sharded tasks.
    """
    if not partitions or partitions == IOfN(0, 1):
        return
    raise ValueError("Paritioned reads are not currently supported")


_ETERNITY_MS = 2**64 - 1


def to_timestamp_ms(input: OpenTimestamp) -> int:
    """Converts a timestamp input type to millis since the Unix epoch."""
    check_type("tiledb_timestamp", input, (int, datetime.datetime))
    if isinstance(input, int):
        timestamp_ms = input
    else:
        # Manually pull out the milliseconds so that nothing funny happens
        # like 12:00:00.300 turning into 299 ms.
        milli_part = input.microsecond // 1000
        input = input.replace(microsecond=0)
        seconds_part = int(input.timestamp())
        timestamp_ms = seconds_part * 1000 + milli_part

    if not 0 <= timestamp_ms <= _ETERNITY_MS:
        raise ValueError("open timestamp must be between 0 (Unix epoch) and 2**64-1 ms")
    return timestamp_ms


def ms_to_datetime(millis: int) -> datetime.datetime:
    """Returns the millisecond timestamp as a timezone-aware UTC datetime.

    This may raise an exception, since millis may be outside the representable
    range for a Python datetime.
    """
    secs, millis = divmod(millis, 1000)
    dt = datetime.datetime.fromtimestamp(secs, tz=datetime.timezone.utc)
    return dt.replace(microsecond=millis * 1000)


def tiledb_timestamp_to_ms(tiledb_timestamp: OpenTimestamp | None) -> int:
    if isinstance(tiledb_timestamp, datetime.datetime):
        tiledb_timestamp = to_timestamp_ms(tiledb_timestamp)
    if tiledb_timestamp is not None and tiledb_timestamp != 0:
        return to_timestamp_ms(tiledb_timestamp)
    return int(time.time() * 1000)


def to_clib_result_order(result_order: ResultOrderStr) -> clib.ResultOrder:
    result_order = ResultOrder(result_order)
    to_clib_result_order = {
        ResultOrder.AUTO: clib.ResultOrder.automatic,
        ResultOrder.ROW_MAJOR: clib.ResultOrder.rowmajor,
        ResultOrder.COLUMN_MAJOR: clib.ResultOrder.colmajor,
    }
    try:
        return to_clib_result_order[result_order]
    except KeyError as ke:
        raise ValueError(f"Invalid result_order: {result_order}") from ke


def from_clib_result_order(result_order: clib.ResultOrder) -> ResultOrder:
    from_clib_result_order = {
        clib.ResultOrder.automatic: ResultOrder.AUTO,
        clib.ResultOrder.rowmajor: ResultOrder.ROW_MAJOR,
        clib.ResultOrder.colmajor: ResultOrder.COLUMN_MAJOR,
    }
    try:
        return from_clib_result_order[result_order]
    except KeyError as ke:
        raise ValueError(f"Invalid clib result_order: {result_order}") from ke


def pa_types_is_string_or_bytes(dtype: pa.DataType) -> bool:
    return bool(
        pa.types.is_large_string(dtype)
        or pa.types.is_large_binary(dtype)
        or pa.types.is_string(dtype)
        or pa.types.is_binary(dtype),
    )


def _cast_domainish(domainish: list[Any]) -> tuple[tuple[object, object], ...]:
    result = []
    for slot in domainish:
        arrow_type = slot[0].type
        if pa.types.is_timestamp(arrow_type):
            result.append(tuple(pa.scalar(to_unix_ts(e), type=arrow_type) for e in slot))
        else:
            result.append(tuple(e.as_py() for e in slot))

    return tuple(result)


def _resolve_futures(unresolved: dict[str, Any], deep: bool = False) -> dict[str, Any]:
    """Resolves any futures found in the dict."""
    resolved = {}
    for k, v in unresolved.items():
        if isinstance(v, Future):
            v = v.result()

        if deep and isinstance(v, dict):
            v = _resolve_futures(v, deep=deep)

        resolved[k] = v

    return resolved


def to_unix_ts(dt: int | pa.TimestampScalar | np.datetime64) -> int:
    if isinstance(dt, pa.TimestampScalar):
        return int(dt.value)
    if isinstance(dt, np.datetime64):
        return int(dt.astype("int64"))
    return dt


class Sentinel:
    """This is used to help detect when a kwarg is supplied or not.  Often,
    kwarg ``foo = None`` is adequate, and should be used. This helps for cases
    when ``None`` is actually a valid option, and we want to distinguish between
    the user passing ``foo=None`` and the user not passing any ``foo`` at all.
    """


MISSING = Sentinel()


def sanitize_key(key: str, data_protocol: DataProtocol) -> str:
    # Encode everything outside of the safe characters set

    if data_protocol == "tiledbv3":
        # Carrara data model supports anything exclusive of '/'
        if "/" in key:
            raise ValueError(f"{key} is not a supported name - must not contain slash (/)")
        sanitized_name = key
    elif data_protocol == "tiledbv2":
        safe_puncuation = "-_.()^!@+={}~'"
        safe_character_set = f"{digits}{ascii_lowercase}{ascii_uppercase}{safe_puncuation}"
        sanitized_name = urllib.parse.quote(key, safe=safe_character_set)
    else:
        raise ValueError(f"Unknown data protocol {data_protocol}")

    # Ensure that the final key is valid
    if sanitized_name in ["..", "."]:
        raise ValueError(f"{key} is not a supported name")

    return sanitized_name


def _df_set_index(
    df: pd.DataFrame,
    default_index_name: str | None = None,
    fallback_index_name: str | None = None,
) -> None:
    if default_index_name is not None:
        # One or both of the following was true:
        # - Original DataFrame had an index name (other than "index") ⇒ that name was written as `OriginalIndexMetadata`
        # - `default_index_name` was provided (e.g. `{obs,var}_id_name` args to `to_anndata`)
        #
        # ⇒ Verify a column with that name exists, and set it as index (keeping its name).
        if default_index_name not in df:
            raise ValueError(f"Requested ID column name {default_index_name} not found in input: {df.keys()}")
        df.set_index(default_index_name, inplace=True)

    else:
        # The assumption here is that the original index was unnamed, and was given a "fallback name" (e.g. "obs_id",
        # "var_id") during ingest that matches the `fallback_index_name` arg here. In this case, we restore that column
        # as index, and remove the name.
        #
        # NOTE: several edge cases result in the outgested DF not matching the original DF; see
        # https://github.com/single-cell-data/TileDB-SOMA/issues/2829.
        if fallback_index_name is not None and fallback_index_name in df:
            df.set_index(fallback_index_name, inplace=True)
            df.index.name = None


def _cast_record_batch(batch: pa.RecordBatch, target_schema: pa.Schema, safe: bool = True) -> pa.RecordBatch:
    """Cast a pyarrow ``RecordBatch`` to another schema.

    ``RecordBatch.cast`` is added in pyarrow==1.16.0. If/when we upgrade our minimum version we can switch
    to the arrow method.

    This method is copied directly from pyarrow.
    """
    if batch.schema.names != target_schema.names:
        raise ValueError(
            f"Target schema's field names are not matching "
            f"the record batch's field names: {batch.schema.names!r}, {target_schema.names!r}"
        )

    newcols = []
    for index in range(batch.num_columns):
        column = batch.column(index)
        field = target_schema.field(index)
        if not field.nullable and column.null_count > 0:
            raise ValueError(f"Casting field {field.name!r} with null values to non-nullable")
        casted = column.cast(field.type, safe=safe)
        newcols.append(casted)
    return pa.RecordBatch.from_arrays(newcols, schema=target_schema)
