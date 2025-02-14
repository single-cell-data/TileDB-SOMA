## Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
##
## Licensed under the MIT License.
#
from __future__ import annotations

from typing import Dict, Union

import pandas as pd
import pyarrow as pa

from tiledbsoma.io.conversions import _prepare_df_for_ingest, df_to_arrow_table

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
    Converts an Arrow schema to a string/string dict, which is easier on the eyes,
    easier to convert from/to JSON for distributed logging, and easier to do del-key on.
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


def _string_dict_from_pandas_dataframe(
    df: pd.DataFrame,
    default_index_name: str,
) -> Dict[str, str]:
    """
    Here we provide compatibility with the ingestor.

    SOMA experiments are indexed by int64 soma_joinid and this is SOMA-only.

    AnnData inputs have a column offered as the index. This can be: named explicitly "obs_id",
    "var_id", etc.; unnamed: adata.obs.index.name is None; named "index".

    In the latter two cases the ingestor allows a rename to the user's choice
    such as "obs_id" and "var_id". Here in the appender pre-check logic, we
    allow the same.
    """

    df = df.head(1).copy()  # since reset_index can be expensive on full data
    _prepare_df_for_ingest(df, default_index_name)
    arrow_table = df_to_arrow_table(df)
    arrow_schema = arrow_table.schema.remove_metadata()
    return _string_dict_from_arrow_schema(arrow_schema)


# Metadata indicating a SOMA DataFrame's original index column name, serialized as a JSON string or `null`.
# SOMA DataFrames are always given a `soma_joinid` index, but we want to be able to outgest a `pd.DataFrame` that is
# identical to the one we ingested, so we store an "original index name" in the DataFrame's metadata.
OriginalIndexMetadata = Union[None, str]
