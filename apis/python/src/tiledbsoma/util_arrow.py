from typing import Iterator, Optional, Union

import numpy as np
import pyarrow as pa
import tiledb

"""
Conversion to/from Arrow and TileDB type systems. Must be capable
of representing full type semantics, and correctly performing a
round trip conversion (eg, T == to_arrow(to_tiledb(T)))

Most primitive types are simple - eg, uint8. Of particular challenge
are datetime/timestamps as TileDB has no distinction between a "datetime" and
a "timedelta". The best Arrow match is TimestampType, as long as that
TimestampType instance does NOT have a timezone set.

Because of our round-trip requirement, all other Arrow temporal types
are unsupported (even though they are just int64 under the covers).
"""
ARROW_TO_TDB = {
    # Dict of types unsupported by to_pandas_dtype, which require overrides.
    # If the value is an instance of Exception, it will be raised.
    #
    # IMPORTANT: ALL non-primitive types supported by TileDB must be in this table.
    #
    pa.string(): np.dtype(
        "S"
    ),  # XXX TODO: temporary work-around until UTF8 support is native. GH #338.
    pa.binary(): np.dtype("S"),
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


def tiledb_type_from_arrow_type(t: pa.DataType) -> Union[type, np.dtype]:
    """
    Given an Arrow type, return the corresponding TileDB type as a Numpy dtype.
    Building block for Arrow-to-TileDB schema translation.

    If type is unsupported, with raise a TypeError exception.

    Parameters
    ----------
    t : pyarrow.DataType
        Arrow DataType instance, eg, pyarrow.int8()

    Returns
    -------
    numpy.dtype
        The numpy dtype corresponding to the ``t`` parameter. ``TypeError`` will
        be raised for unsupported types.
    """
    if t in ARROW_TO_TDB:
        arrow_type = ARROW_TO_TDB[t]
        if isinstance(arrow_type, Exception):
            raise arrow_type
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


def get_arrow_type_from_tiledb_dtype(tiledb_dtype: np.dtype) -> pa.DataType:
    """
    TODO: COMMENT
    """
    if tiledb_dtype.name == "bytes":
        # XXX TODO: temporary work-around until UTF8 support is native. GH #338.
        return pa.string()
    else:
        return pa.from_numpy_dtype(tiledb_dtype)


def get_arrow_schema_from_tiledb_uri(
    tiledb_uri: str, ctx: Optional[tiledb.Ctx] = None
) -> pa.Schema:
    """
    TODO: COMMENT
    """
    with tiledb.open(tiledb_uri, ctx=ctx) as A:
        arrow_schema_dict = {}

        dom = A.schema.domain
        for i in range(dom.ndim):
            dim = dom.dim(i)
            name = dim.name
            if name == "":
                name = "unnamed"
            arrow_schema_dict[name] = get_arrow_type_from_tiledb_dtype(dim.dtype)

        for i in range(A.schema.nattr):
            attr = A.schema.attr(i)
            name = attr.name
            if name == "":
                name = "unnamed"
            arrow_schema_dict[name] = get_arrow_type_from_tiledb_dtype(attr.dtype)

    return pa.schema(arrow_schema_dict)


def ascii_to_unicode_pyarrow_readback(record_batch: pa.RecordBatch) -> pa.RecordBatch:
    """
    Implements the 'decode on read' part of our ASCII/Unicode logic
    """
    # TODO: COMMENT/LINK HEAVILY
    names = [ofield.name for ofield in record_batch.schema]
    new_fields = []
    for name in names:
        old_field = record_batch[name]
        if isinstance(old_field, pa.LargeBinaryArray):
            nfield = pa.array(
                [element.as_py().decode("utf-8") for element in old_field]
            )
            new_fields.append(nfield)
        else:
            new_fields.append(old_field)
    return pa.RecordBatch.from_arrays(new_fields, names=names)


def concat_batches(batch_generator: Iterator[pa.RecordBatch]) -> pa.RecordBatch:
    """
    Iterates a generator of ``pyarrow.RecordBatch`` (e.g. ``SOMADataFrame.read``) and returns a concatenation of all the record batches found. The nominal use is to simply unit-test cases.
    """
    batches = []
    for batch in batch_generator:
        batches.append(batch)
    assert len(batches) > 0
    names = [field.name for field in batches[0].schema]
    arrays = []
    for name in names:
        array = pa.concat_arrays([batch[name] for batch in batches])
        arrays.append(array)
    return pa.RecordBatch.from_arrays(arrays, names=names)
