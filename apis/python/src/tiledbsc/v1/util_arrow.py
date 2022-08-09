from typing import Generator, Union

import numpy as np
import pyarrow as pa


def tiledb_type_from_arrow_type(t: pa.DataType) -> Union[type, np.dtype]:
    """
    Building block for Arrow-to-TileDB schema translation.
    """
    if t == pa.string():
        # pyarrow's to_pandas_dtype maps pa.string() to dtype object which
        # isn't acceptable to tiledb -- we must say str.
        # XXX COMMENT return str
        return np.dtype("S")
    else:
        # mypy says:
        # Returning Any from function declared to return "type"  [no-any-return]
        return t.to_pandas_dtype()  # type: ignore


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


def concat_batches(batch_generator: Generator) -> pa.RecordBatch:
    """
    Iterates a generator of `pyarrow.RecordBatch` (e.g. `SOMADataFrame.read`) and returns a
    concatenation of all the record batches found. The nominal use is to simply unit-test cases.
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
