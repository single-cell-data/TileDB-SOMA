# Copyright (c) 2021-2023 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2023 TileDB, Inc.
#
# Licensed under the MIT License.

import datetime
import json
import pathlib
import time
import urllib.parse
from itertools import zip_longest
from typing import Any, Dict, List, Mapping, Optional, Tuple, Type, TypeVar, Union, cast

import pyarrow as pa
import somacore
from somacore import options

from . import pytiledbsoma as clib
from ._types import OpenTimestamp, Slice, is_slice_of
from .options._tiledb_create_write_options import (
    TileDBCreateOptions,
    _ColumnConfig,
    _DictFilterSpec,
)

_JSONFilter = Union[str, Dict[str, Union[str, Union[int, float]]]]
_JSONFilterList = Union[str, List[_JSONFilter]]


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
    return "%s TIME %.3f seconds" % (message, time.time() - start_stamp)


def is_local_path(path: str) -> bool:
    if path.startswith("file://"):
        return True
    if "://" in path:
        return False
    return True


def make_relative_path(uri: str, relative_to: str) -> str:
    """Returns a URI relative to another URI. If not possible, raise a ValueError.

    This function assumes that the URI scheme follows posix path conventions
    and only contains a scheme, netloc and path. It does not handle query params,
    fragments, etc.

    The default scheme, if one is not specified, is assumed to be `file`.
    """
    p_uri = urllib.parse.urlparse(uri)
    p_relative_to = urllib.parse.urlparse(relative_to)

    uri_scheme = p_uri.scheme if p_uri.scheme != "" else "file"
    relative_to_scheme = p_relative_to.scheme if p_relative_to.scheme != "" else "file"
    if uri_scheme != relative_to_scheme:
        raise ValueError(
            "Unable to make relative path between URIs with different scheme"
        )

    relpath = pathlib.PurePath(p_uri.path).relative_to(p_relative_to.path).as_posix()
    return relpath


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

    if p_base.scheme == "" or p_base.scheme == "file":
        # if a file path, just use pathlib.
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
    if slc.stop < slc.start:
        raise ValueError(
            f"slice start ({slc.start!r}) must be <= slice stop ({slc.stop!r})"
        )


_T = TypeVar("_T")


class NonNumericDimensionError(TypeError):
    """Raised when trying to get a numeric range for a non-numeric dimension."""


def slice_to_numeric_range(
    slc: Slice[Any], domain: Tuple[_T, _T]
) -> Optional[Tuple[_T, _T]]:
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

    if stop < start:
        # With the above, we have guaranteed that at least one bound will
        # include the domain.  If we get here, that means that the other bound
        # never included it (e.g. stop == slc.stop < domain_start == start).
        raise ValueError(
            f"slice [{slc.start!r}:{slc.stop!r}] does not overlap"
            f" [{domain_start!r}:{domain_stop!r}]"
        )

    return start, stop


def dense_indices_to_shape(
    coords: options.DenseNDCoords,
    array_shape: Tuple[int, ...],
    result_order: somacore.ResultOrder,
) -> Tuple[int, ...]:
    """Given a subarray index specified as a tuple of per-dimension slices or scalars
    (e.g., ``([:], 1, [1:2])``), and the shape of the array, return the shape of
    the subarray. Note that the number of coordinates may be less than or equal
    to the number of dimensions in the array.
    """
    if len(coords) > len(array_shape):
        raise ValueError(
            f"coordinate length ({len(coords)}) must be <="
            f" array dimension count ({len(array_shape)})"
        )

    shape = tuple(
        dense_index_to_shape(coord, extent)
        for coord, extent in zip_longest(coords, array_shape)
    )
    if result_order is somacore.ResultOrder.ROW_MAJOR:
        return shape
    return tuple(reversed(shape))


def dense_index_to_shape(coord: options.DenseCoord, array_length: int) -> int:
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
        if coord.stop is None:
            stop = array_length
        else:
            stop = min(coord.stop + 1, array_length)
        return stop - (coord.start or 0)

    raise TypeError(f"coordinate {coord} must be integer or integer slice")


def check_type(
    name: str,
    actual_value: Any,
    expected_types: Tuple[Type[Any], ...],
) -> None:
    """Verifies the type of an argument, or produces a useful error message."""
    if not isinstance(actual_value, expected_types):
        if len(expected_types) == 1:
            raise TypeError(
                f"expected {name} argument to be of type {expected_types[0]}; got {type(actual_value)}"
            )
        raise TypeError(
            f"expected {name} argument to be one of {expected_types!r}; got {type(actual_value)}"
        )


def check_unpartitioned(partitions: Optional[options.ReadPartitions]) -> None:
    """Ensures that we're not being asked for a partitioned read.

    Because we currently don't support partitioned reads, we should reject all
    reads that request partitions to avoid giving the user duplicate data across
    sharded tasks.
    """
    if not partitions or partitions == options.IOfN(0, 1):
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


def to_clib_result_order(result_order: options.ResultOrderStr) -> clib.ResultOrder:
    result_order = options.ResultOrder(result_order)
    to_clib_result_order = {
        options.ResultOrder.AUTO: clib.ResultOrder.automatic,
        options.ResultOrder.ROW_MAJOR: clib.ResultOrder.rowmajor,
        options.ResultOrder.COLUMN_MAJOR: clib.ResultOrder.colmajor,
    }
    try:
        return to_clib_result_order[result_order]
    except KeyError as ke:
        raise ValueError(f"Invalid result_order: {result_order}") from ke


def pa_types_is_string_or_bytes(dtype: pa.DataType) -> bool:
    return bool(
        pa.types.is_large_string(dtype)
        or pa.types.is_large_binary(dtype)
        or pa.types.is_string(dtype)
        or pa.types.is_binary(dtype)
    )


def cast_values_to_target_schema(values: pa.Table, schema: pa.Schema) -> pa.Table:
    """
    When writing data to a SOMAArray, the values that the user passes in may not
    match the schema on disk. Cast the values to the correct dtypes.
    """
    # Ensure fields are in the correct order
    # target_schema = []
    # for input_field in values.schema:
    #     target_schema.append(schema.field(input_field.name))

    # return values.cast(pa.schema(target_schema, values.schema.metadata))
    return values


def build_clib_platform_config(
    platform_config: Optional[options.PlatformConfig],
) -> clib.PlatformConfig:
    """
    Copy over Python PlatformConfig values to the C++ clib.PlatformConfig
    """
    plt_cfg = clib.PlatformConfig()

    if platform_config is None:
        return plt_cfg

    ops = TileDBCreateOptions.from_platform_config(platform_config)
    plt_cfg.dataframe_dim_zstd_level = ops.dataframe_dim_zstd_level
    plt_cfg.sparse_nd_array_dim_zstd_level = ops.sparse_nd_array_dim_zstd_level
    plt_cfg.dense_nd_array_dim_zstd_level = ops.dense_nd_array_dim_zstd_level
    plt_cfg.write_X_chunked = ops.write_X_chunked
    plt_cfg.goal_chunk_nnz = ops.goal_chunk_nnz
    plt_cfg.capacity = ops.capacity
    plt_cfg.offsets_filters = _build_filter_list(ops.offsets_filters)
    plt_cfg.validity_filters = _build_filter_list(ops.validity_filters)
    plt_cfg.allows_duplicates = ops.allows_duplicates
    plt_cfg.tile_order = ops.tile_order
    plt_cfg.cell_order = ops.cell_order
    plt_cfg.dims = _build_column_config(ops.dims)
    plt_cfg.attrs = _build_column_config(ops.attrs)
    return plt_cfg


def _build_column_config(col: Optional[Mapping[str, _ColumnConfig]]) -> str:
    column_config: Dict[str, Dict[str, Union[_JSONFilterList, int]]] = dict()

    if col is None:
        return ""

    for k in col:
        dikt: Dict[str, Union[_JSONFilterList, int]] = {}
        if col[k].filters is not None:
            dikt["filters"] = _build_filter_list(col[k].filters, False)
        if col[k].tile is not None:
            dikt["tile"] = cast(int, col[k].tile)
        if len(dikt) != 0:
            column_config[k] = dikt
    return json.dumps(column_config)


def _build_filter_list(
    filters: Optional[Tuple[_DictFilterSpec, ...]], return_json: bool = True
) -> _JSONFilterList:
    _convert_filter = {
        "GzipFilter": "GZIP",
        "ZstdFilter": "ZSTD",
        "LZ4Filter": "LZ4",
        "Bzip2Filter": "BZIP2",
        "RleFilter": "RLE",
        "DeltaFilter": "DELTA",
        "DoubleDeltaFilter": "DOUBLE_DELTA",
        "BitWidthReductionFilter": "BIT_WIDTH_REDUCTION",
        "BitShuffleFilter": "BITSHUFFLE",
        "ByteShuffleFilter": "BYTESHUFFLE",
        "PositiveDeltaFilter": "POSITIVE_DELTA",
        "ChecksumMD5Filter": "CHECKSUM_MD5",
        "ChecksumSHA256Filter": "CHECKSUM_SHA256",
        "DictionaryFilter": "DICTIONARY_ENCODING",
        "FloatScaleFilter": "SCALE_FLOAT",
        "XORFilter": "XOR",
        "WebpFilter": "WEBP",
        "NoOpFilter": "NOOP",
    }

    _convert_option = {
        "GZIP": {"level": "COMPRESSION_LEVEL"},
        "ZSTD": {"level": "COMPRESSION_LEVEL"},
        "LZ4": {"level": "COMPRESSION_LEVEL"},
        "BZIP2": {"level": "COMPRESSION_LEVEL"},
        "RLE": {"level": "COMPRESSION_LEVEL"},
        "DELTA": {
            "level": "COMPRESSION_LEVEL",
            "reinterp_dtype": "COMPRESSION_REINTERPRET_DATATYPE",
        },
        "DOUBLE_DELTA": {
            "level": "COMPRESSION_LEVEL",
            "reinterp_dtype": "COMPRESSION_REINTERPRET_DATATYPE",
        },
        "DICTIONARY_ENCODING": {"level": "COMPRESSION_LEVEL"},
        "BIT_WIDTH_REDUCTION": {"window": "BIT_WIDTH_MAX_WINDOW"},
        "POSITIVE_DELTA": {"window": "POSITIVE_DELTA_MAX_WINDOW"},
        "SCALE_FLOAT": {
            "factor": "SCALE_FLOAT_FACTOR",
            "offset": "SCALE_FLOAT_OFFSET",
            "bytewidth": "SCALE_FLOAT_BYTEWIDTH",
        },
        "WEBP": {
            "input_format": "WEBP_INPUT_FORMAT",
            "quality": "WEBP_QUALITY",
            "lossless": "WEBP_LOSSLESS",
        },
    }

    if filters is None:
        return ""

    filter: _JSONFilter
    filter_list: List[_JSONFilter] = []

    for info in filters:
        if len(info) == 1:
            filter = _convert_filter[cast(str, info["_type"])]
        else:
            filter = dict()
            for option_name, option_value in info.items():
                filter_name = _convert_filter[cast(str, info["_type"])]
                if option_name == "_type":
                    filter["name"] = filter_name
                else:
                    filter[_convert_option[filter_name][option_name]] = cast(
                        Union[float, int], option_value
                    )
        filter_list.append(filter)
    return json.dumps(filter_list) if return_json else filter_list
