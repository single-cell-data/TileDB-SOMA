# This ABSOLUTELY MUST be imported before pyarrow to avoid abort-traps on MacOS
# with any pyarrow >= 13.  And the pre-commit hook, isort, et al. MUST NOT be
# allowed to reorder this.  See also
# https://github.com/single-cell-data/TileDB-SOMA/issues/1926
# https://github.com/single-cell-data/TileDB-SOMA/issues/1849
import tiledb  # noqa E402

import pyarrow as pa
from typeguard import install_import_hook

install_import_hook("tiledbsoma")

"""Types supported in a SOMA*NDArray """
NDARRAY_ARROW_TYPES_SUPPORTED = [
    pa.bool_(),
    pa.int8(),
    pa.int16(),
    pa.int32(),
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
]

"""Primitive types NOT supported in a SOMA*NDArray """
NDARRAY_ARROW_TYPES_NOT_SUPPORTED = [
    pa.float16(),
    pa.date32(),
    pa.time32("s"),
    pa.duration("s"),
    pa.date64(),
    pa.time64("us"),
    pa.time64("ns"),
    pa.binary(),
    pa.binary(10),
    pa.large_binary(),
    pa.large_string(),
    pa.decimal128(1),
    pa.decimal128(38),
    pa.list_(pa.int8()),
    pa.large_list(pa.bool_()),
    pa.map_(pa.string(), pa.int32()),
    pa.struct([("f1", pa.int32()), ("f2", pa.string())]),
    pa.dictionary(pa.int32(), pa.string()),
]
