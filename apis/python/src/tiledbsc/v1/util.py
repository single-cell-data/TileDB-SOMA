import re

import pyarrow as pa

SOMA_OBJECT_TYPE_METADATA_KEY = "soma_object_type"
SOMA_ENCODING_VERSION_METADATA_KEY = "soma_encoding_version"
SOMA_ENCODING_VERSION = "1"


def tiledb_type_from_arrow_type(t: pa.DataType) -> type:
    """
    TODO
    """
    if t == pa.string():
        # pyarrow's to_pandas_dtype maps pa.string() to dtype object which
        # isn't acceptable to tiledb -- we must say str.
        return str
    else:
        # mypy says:
        # Returning Any from function declared to return "type"  [no-any-return]
        return t.to_pandas_dtype()  # type: ignore


def is_tiledb_creation_uri(uri: str) -> bool:
    return bool(re.match("^tiledb://.*s3://.*$", uri))
