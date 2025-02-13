## Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
##
## Licensed under the MIT License.
#
from __future__ import annotations

from typing import Dict, Union

import pandas as pd
import pyarrow as pa

from tiledbsoma.io.conversions import df_to_arrow

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
    arrow_table = df_to_arrow(df)
    arrow_schema = arrow_table.schema.remove_metadata()
    return _string_dict_from_arrow_schema(arrow_schema)


# Metadata indicating a SOMA DataFrame's original index column name, serialized as a JSON string or `null`.
# SOMA DataFrames are always given a `soma_joinid` index, but we want to be able to outgest a `pd.DataFrame` that is
# identical to the one we ingested, so we store an "original index name" in the DataFrame's metadata.
OriginalIndexMetadata = Union[None, str]


def _prepare_df_for_ingest(
    df: pd.DataFrame, id_column_name: str | None
) -> OriginalIndexMetadata:
    """Prepare a `pd.DataFrame` for persisting as a SOMA DataFrame: demote its index to a column (to make way for a
    required `soma_joinid` index), and compute and return metadata for restoring the index column and name later (on
    outgest).

    If `df.index` has a name (and it's not "index", which is taken to be a default/unset value):
    - `df.index.name` takes precedence over the `id_column_name` arg: the index will be reset to an eponymous column.
    - That original `df.index.name` will be logged as `OriginalIndexMetadata` (for promotion back to index on outgest).

    In this case, the overall round trip is basically just:
    - `reset_index` on ingest (demote index to eponymous column).
    - `set_index` on outgest (restore column to index, with its original name).

    Otherwise (index name is `None` or "index"):
    - A fallback name (`id_column_name` if provided, "index" otherwise) is used for the column that the index becomes.
    - The returned `OriginalIndexMetadata` will be `None`.

    There are several edge cases (detailed below and in `test_dataframe_io_roundtrips.py` and
    https://github.com/single-cell-data/TileDB-SOMA/issues/2829) where the index, its name, or a specific column are not
    restored properly on outgest. For now, all such behavior is preserved, for backwards compatibility, but we should
    look into ways of improving these "round-trip mutation" cases. See
    https://github.com/single-cell-data/TileDB-SOMA/issues/2829 for more info.
    """
    use_existing_index = df.index.name is not None and df.index.name != "index"

    original_index_name = None
    if use_existing_index:
        original_index_name = df.index.name

    df.reset_index(inplace=True)
    if id_column_name is not None:
        if id_column_name in df:
            if "index" in df:
                # The assumption here is that the column named "index" was previously an unnamed `df.index`, and
                # `id_column_name` was already a column (per the grandparent `if` above). In this case, we drop the
                # original unnamed `df.index`.
                # TODO: This prevents outgesting the same DataFrame we ingested. We should fix it; see
                #  https://github.com/single-cell-data/TileDB-SOMA/issues/2829.
                #
                # Also note: if the DataFrame already had columns named "index" and `id_column_name`, the original
                # `df.index` will have been "reset" to a column named `level_0`, and we end up just dropping the column
                # named "index" here.
                #
                # Another version of this occurs when the original DataFrame has `df.index.name == id_column_name` and a
                # column named "index". In this case, the index will have been "reset" to a column named
                # `id_column_name` above, which then satisfies the grendparent `if`'s predicate, and causes us to drop
                # the column named "index" here.
                df.drop(columns=["index"], inplace=True)
        else:
            # If `id_column_name` was passed, and is not already a column in the DataFrame, we assume the original index
            # was "reset" to a column named "index" (by `reset_index` above), and we rename that column to
            # `id_column_name`, so that `id_column_name` matches the name of a column representing the original
            # DataFrame's index.
            #
            # NOTE: the assumption above can break in a few ways:
            # 1. The original DataFrame index has a name other than "index" or `id_column_name`…
            #    a. and there is a column named "index" ⇒ that column will be renamed to `id_column_name`
            #    b. and there is no column named "index" ⇒ the rename below is a no-op (outgest currently restores the
            #       original DataFrame in this case)
            # 2. The original DataFrame has a column named "index":
            #    - That column will become `df.index` on outgest, and acquire the original `df.index.name` as its name.
            #    - The original index will end up as a column, on outgest:
            #      - If it had a name, the column will have that name.
            #      - Otherwise, it will end up as a column named e.g. `level_0` (or `level_1`, if a column named
            #        `level_0` already exists, etc.)
            #
            # See https://github.com/single-cell-data/TileDB-SOMA/issues/2829 for more info.
            df.rename(columns={"index": id_column_name}, inplace=True)

    return original_index_name
