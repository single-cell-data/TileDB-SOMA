import numpy as np
import pyarrow as pa
import pytest

from tiledbsoma.util_arrow import (
    get_arrow_type_from_tiledb_dtype,
    tiledb_type_from_arrow_type,
)

"""Arrow types we expect to work"""
SUPPORTED_ARROW_TYPES = [
    pa.bool_(),
    pa.int8(),
    pa.int16(),
    pa.int32(),
    pa.int16(),
    pa.uint8(),
    pa.uint16(),
    pa.uint32(),
    pa.uint64(),
    pa.float32(),
    pa.float64(),
    pa.timestamp("s"),
    pa.timestamp("ms"),
    pa.timestamp("us"),
    pa.timestamp("ns"),
    # We use Arrow's large_string for ASCII and, ultimately, for Unicode as well
    # https://github.com/single-cell-data/TileDB-SOMA/issues/99
    # https://github.com/single-cell-data/TileDB-SOMA/pull/359
    # https://github.com/single-cell-data/TileDB-SOMA/issues/274
    pa.large_string(),
    pa.large_binary(),
]

"""Arrow types we expect to auto-promote"""
PROMOTED_ARROW_TYPES = [
    (pa.string(), pa.large_string()),
    # XXX (pa.binary(), pa.large_binary()),
]


"""Arrow types we expect to fail"""
UNSUPPORTED_ARROW_TYPES = [
    pa.null(),
    pa.date64(),
    pa.time64("us"),
    pa.time64("ns"),
    pa.float16(),
    pa.date32(),
    pa.time32("s"),
    pa.time32("ms"),
    pa.duration("s"),
    pa.duration("ms"),
    pa.duration("us"),
    pa.duration("ns"),
    pa.month_day_nano_interval(),
    # We use Arrow's large_string for ASCII and, ultimately, for Unicode as well
    # https://github.com/single-cell-data/TileDB-SOMA/issues/99
    # https://github.com/single-cell-data/TileDB-SOMA/pull/359
    # https://github.com/single-cell-data/TileDB-SOMA/issues/274
    pa.string(),
    pa.binary(),
    pa.binary(10),
    pa.decimal128(1),
    pa.decimal128(38),
    pa.list_(pa.int8()),
    pa.large_list(pa.bool_()),
    pa.map_(pa.string(), pa.int32()),
    pa.struct([("f1", pa.int32()), ("f2", pa.string())]),
    pa.dictionary(pa.int32(), pa.string()),
]


@pytest.mark.parametrize("arrow_type", SUPPORTED_ARROW_TYPES)
def test_arrow_types_supported(arrow_type):
    """Verify round-trip conversion of types"""
    # if pa.types.is_binary(arrow_type):
    # pytest.xfail("Awaiting UTF-8 support - see issue #274")

    tdb_dtype = tiledb_type_from_arrow_type(arrow_type)
    assert (
        isinstance(tdb_dtype, np.dtype) or tdb_dtype == "ascii" or tdb_dtype == "bytes"
    )
    arrow_rt_type = get_arrow_type_from_tiledb_dtype(tdb_dtype)
    assert isinstance(arrow_rt_type, pa.DataType)
    assert arrow_type == arrow_rt_type


@pytest.mark.parametrize("arrow_from_to_pair", PROMOTED_ARROW_TYPES)
def test_arrow_types_promoted(arrow_from_to_pair):
    """Verify round-trip conversion of types"""
    arrow_from_type = arrow_from_to_pair[0]
    arrow_to_type = arrow_from_to_pair[1]

    tdb_dtype = tiledb_type_from_arrow_type(arrow_from_type)
    assert (
        isinstance(tdb_dtype, np.dtype) or tdb_dtype == "ascii" or tdb_dtype == "bytes"
    )
    arrow_rt_type = get_arrow_type_from_tiledb_dtype(tdb_dtype)
    assert isinstance(arrow_rt_type, pa.DataType)
    assert arrow_to_type == arrow_rt_type


@pytest.mark.parametrize("arrow_type", UNSUPPORTED_ARROW_TYPES)
def test_arrow_types_unsupported(arrow_type):
    """Verify correct error for unsupported types"""
    with pytest.raises(TypeError):
        tiledb_type_from_arrow_type(arrow_type, match=r".*unsupported type.*")
