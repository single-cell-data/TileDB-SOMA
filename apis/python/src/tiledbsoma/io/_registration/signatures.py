## Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
##
## Licensed under the MIT License.
#
from __future__ import annotations

from typing import Dict, Union

import pyarrow as pa

_EQUIVALENCES = {
    "large_string": "string",
    "large_binary": "binary",
}


def _stringify_type(t: pa.DataType) -> str:
    """
    Turns an Arrow data type into a stringi more suitable for logging error messages to users in a
    distributed-computing/distributed-logging environment.

    As noted in the Signature class, we pre-check logic from the ingestor.  As detailed elsewhere,
    Arrow string and large_string must map to TileDB string, which is large-only. Thus string and
    large_string form an equivalence class. Similarly for Arrow binary and large_binary.
    """
    str_t = str(t)
    return _EQUIVALENCES.get(str_t, str_t)


def _string_dict_from_arrow_schema(schema: pa.Schema) -> Dict[str, str]:
    """
    Converts an Arrow schema to a string/string dict.

    This is easier on the eyes, easier to convert from/to JSON for distributed logging,
    and easier to do del-key on.
    """
    retval = {}
    for name in schema.names:
        arrow_type = schema.field(name).type
        if pa.types.is_dictionary(arrow_type):
            arrow_type = arrow_type.index_type
        retval[name] = _stringify_type(arrow_type)
    # The soma_joinid field is specific to SOMA data but does not exist in AnnData/H5AD.  When we
    # pre-check an AnnData/H5AD input to see if it's appendable to an existing SOMA experiment, we
    # must not punish the AnnData/H5AD input for it not having a soma_joinid column in its obs and
    # var.
    retval.pop("soma_joinid", None)
    return retval


# Metadata indicating a SOMA DataFrame's original index column name, serialized as a JSON string or `null`.
# SOMA DataFrames are always given a `soma_joinid` index, but we want to be able to outgest a `pd.DataFrame` that is
# identical to the one we ingested, so we store an "original index name" in the DataFrame's metadata.
OriginalIndexMetadata = Union[None, str]
