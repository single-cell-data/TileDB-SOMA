"""
Conversion to/from Arrow and TileDB type systems. Must be capable
of representing full type semantics, and correctly performing a
round trip conversion (e.g., T == to_arrow(to_tiledb(T)))

Most primitive types are simple -- e.g., uint8. Of particular challenge
are datetime/timestamps as TileDB has no distinction between a "datetime" and
a "timedelta". The best Arrow match is TimestampType, as long as that
TimestampType instance does NOT have a timezone set.

Because of our round-trip requirement, all other Arrow temporal types
are unsupported (even though they are just int64 under the covers).

We auto-promote Arrow's string and binary to large_string and large_binary,
respectively, as this is what TileDB stores -- a sequence of bytes preceded
by a 64-bit (not 32-bit) length int.

DataFrame-specific note: currently (as of 2.14), TileDB does not support
Unicode array dimensions. All Arrow string types used in a DataFrame index
columns (i.e., TileDB dimension) are coerced to ASCII. This equirement for
ASCII-only dimensions will be relaxed in a future release. Unicode/UTF-8 is
fully supported in SOMA DataFrame non-indexed columns.
"""

from typing import Any, Dict, Union

import numpy as np
import numpy.typing as npt
import pyarrow as pa
import tiledb

_ARROW_TO_TDB_ATTR: Dict[Any, Union[str, TypeError]] = {
    pa.string(): "U1",
    pa.large_string(): "U1",
    pa.binary(): "bytes",
    pa.large_binary(): "bytes",
    pa.timestamp("s"): "datetime64[s]",
    pa.timestamp("ms"): "datetime64[ms]",
    pa.timestamp("us"): "datetime64[us]",
    pa.timestamp("ns"): "datetime64[ns]",
    #
    # Unsupported types in TileDB type system
    pa.float16(): TypeError("float16 - unsupported type (use float32)"),
    pa.date32(): TypeError("32-bit date - unsupported type (use TimestampType)"),
    pa.date64(): TypeError("64-bit date - unsupported type (use TimestampType)"),
}
"""
Dict of types unsupported by to_pandas_dtype, which require overrides for
use in TileDB Attributes (aka DataFrame non-indexe columns).

If the value is an instance of Exception, it will be raised.

IMPORTANT: ALL non-primitive types supported by TileDB must be in this table.
"""

# Same as _ARROW_TO_TDB_ATTR, but used for DataFrame indexed columns, aka TileDB Dimensions.
# Any type system differences from the base-case Attr should be added here.
_ARROW_TO_TDB_DIM: Dict[Any, Union[str, TypeError]] = _ARROW_TO_TDB_ATTR.copy()
"""
Same as _ARROW_TO_TDB_ATTR, but used for DataFrame indexed columns, aka TileDB Dimensions.
Any type system differences from the base-case Attr should be added here.
"""
_ARROW_TO_TDB_DIM.update(
    {
        pa.string(): "ascii",  # TODO: temporary work-around until Dimension UTF8 support is available.
        pa.large_string(): "ascii",  # TODO: temporary work-around until Dimension UTF8 support is available.
    }
)


def tiledb_type_from_arrow_type(
    t: pa.DataType, is_indexed_column: bool = False
) -> npt.DTypeLike:
    """
    Given an Arrow type, return the corresponding TileDB type as a NumPy dtype.
    Building block for Arrow-to-TileDB schema translation.

    TileDB currently has different Unicode handling for dimensions and attributes.
    Set the ``is_dimension`` parameter to True for indexed-column (AKA dimension)
    rules, which currently requires all strings to be ASCII.

    If type is unsupported, with raise a TypeError exception.

    Parameters
    ----------
    t : pyarrow.DataType
        Arrow DataType instance, e.g., pyarrow.int8()
    is_indexed_column : bool
        Use TileDB dimension type conversion rules.

    Returns
    -------
    numpy.dtype
        The numpy dtype corresponding to the ``t`` parameter. ``TypeError`` will
        be raised for unsupported types.
    """
    arrow_to_tdb = _ARROW_TO_TDB_DIM if is_indexed_column else _ARROW_TO_TDB_ATTR
    if t in arrow_to_tdb:
        arrow_type = arrow_to_tdb[t]
        if isinstance(arrow_type, Exception):
            raise arrow_type
        if arrow_type in ["ascii", "bytes"]:
            return arrow_type
        return np.dtype(arrow_type)

    if not pa.types.is_primitive(t):
        raise TypeError(f"Type {str(t)} - unsupported type")
    if pa.types.is_timestamp(t):
        raise TypeError("TimeStampType - unsupported type (timezone not supported)")
    if pa.types.is_time32(t):
        raise TypeError("Time64Type - unsupported type (use TimestampType)")
    if pa.types.is_time64(t):
        raise TypeError("Time32Type - unsupported type (use TimestampType)")
    if pa.types.is_duration(t):
        raise TypeError("DurationType - unsupported type (use TimestampType)")

    # else lets try the default conversion path
    try:
        # Must force into a dtype to catch places where the Pandas type
        # system has extra information that can't be expressed
        return np.dtype(t.to_pandas_dtype())
    except NotImplementedError as exc:
        raise TypeError("Unsupported Arrow type") from exc


def arrow_type_from_tiledb_dtype(
    tiledb_dtype: npt.DTypeLike, bytes_are_ascii: bool = True
) -> pa.DataType:
    """
    Maps a TileDB dtype (``'bytes'``, ``'ascii'``, or an ``np.dtype``) to an Arrow type.  Note that
    when we read tiledb schema off storage, ``ascii`` and ``bytes`` both have ``dtype`` of `"S"`
    which is equal to ``bytes`` -- so, the caller should disambgiuate.
    """
    if tiledb_dtype == "bytes":
        if bytes_are_ascii:
            return pa.large_string()
        else:
            return pa.large_binary()
    elif tiledb_dtype == "ascii" or tiledb_dtype == str:
        return pa.large_string()
    else:
        return pa.from_numpy_dtype(tiledb_dtype)


def tiledb_schema_to_arrow(tdb_schema: tiledb.ArraySchema) -> pa.Schema:
    arrow_schema_dict = {}
    dom = tdb_schema.domain
    for i in range(dom.ndim):
        dim = dom.dim(i)
        name = dim.name
        if name == "":
            name = "unnamed"
        arrow_schema_dict[name] = arrow_type_from_tiledb_dtype(dim.dtype)

    for i in range(tdb_schema.nattr):
        attr = tdb_schema.attr(i)
        name = attr.name
        if name == "":
            name = "unnamed"
        arrow_schema_dict[name] = arrow_type_from_tiledb_dtype(attr.dtype, attr.isascii)

    return pa.schema(arrow_schema_dict)
