# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

from __future__ import annotations

import json
from collections.abc import Mapping
from typing import Union, cast

from tiledbsoma import pytiledbsoma as clib
from tiledbsoma._core_options import PlatformConfig

from ._tiledb_create_write_options import TileDBCreateOptions, _ColumnConfig, _DictFilterSpec

_JSONFilter = Union[str, dict[str, Union[str, Union[int, float]]]]
_JSONFilterList = Union[str, list[_JSONFilter]]


def build_clib_platform_config(
    platform_config: PlatformConfig | None,
) -> clib.PlatformConfig:
    """Copy over Python PlatformConfig values to the C++ clib.PlatformConfig."""
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


def _build_column_config(col: Mapping[str, _ColumnConfig] | None) -> str:
    column_config: dict[str, dict[str, _JSONFilterList | int]] = {}

    if col is None:
        return ""

    for k in col:
        dikt: dict[str, _JSONFilterList | int] = {}
        if col[k].filters is not None:
            dikt["filters"] = _build_filter_list(col[k].filters, False)
        if col[k].tile is not None:
            dikt["tile"] = cast("int", col[k].tile)
        if len(dikt) != 0:
            column_config[k] = dikt
    return json.dumps(column_config)


def _build_filter_list(filters: tuple[_DictFilterSpec, ...] | None, return_json: bool = True) -> _JSONFilterList:
    convert_filter_ = {
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

    convert_option_ = {
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
    filter_list: list[_JSONFilter] = []

    for info in filters:
        if len(info) == 1:
            filter = convert_filter_[cast("str", info["_type"])]
        else:
            filter = {}
            for option_name, option_value in info.items():
                filter_name = convert_filter_[cast("str", info["_type"])]
                if option_name == "_type":
                    filter["name"] = filter_name
                else:
                    filter[convert_option_[filter_name][option_name]] = cast("Union[float, int]", option_value)
        filter_list.append(filter)
    return json.dumps(filter_list) if return_json else filter_list
