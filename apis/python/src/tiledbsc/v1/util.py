import pyarrow as pa


def tiledb_type_from_arrow_type(t: pa.DataType) -> type:
    """
    TODO
    """
    if t == pa.string():
        # pyarrow's to_pandas_dtype maps pa.string() to dtype object which
        # isn't acceptable to tiledb -- we must say str.
        return str
    else:
        return t.to_pandas_dtype()
