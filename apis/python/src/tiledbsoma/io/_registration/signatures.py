import json
from typing import Dict, Optional, Union

import anndata as ad
import attrs
import pandas as pd
import pyarrow as pa
from typing_extensions import Self

import tiledbsoma
import tiledbsoma.logging
from tiledbsoma._arrow_types import df_to_arrow
from tiledbsoma.io._util import read_h5ad  # Allow us to read over S3 in backed mode
from tiledbsoma.options import SOMATileDBContext

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
    df: pd.DataFrame, id_column_name: Optional[str]
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


@attrs.define(kw_only=True)
class Signature:
    """
    This is support for compatibility pre-check for append-mode SOMA ingestion.

    If a SOMA experiment already exists and the user wants to append another AnnData/H5AD to it ---
    or, if no SOMA experiment exists yet but the user has two or more AnnData/H5AD inputs --- we
    provide a fail-fast schema-compatibility check.  Use of this pre-check ensures that we flag
    non-appendable data before any data writes start. In particular, we avoid having a SOMA
    experiment half-appended to.

    At present we require that all schemas are identical, both within the ingestor and this
    pre-checker. This includes ``obs`` and ``var`` field names and exact dtypes, as well as ``X``,
    ``obsm``, and ``varm`` dtypes.

    Later, we can relax constraints: if say the SOMA experiment has ``float64`` dtype for ``X``
    and an H5AD append has ``float32``, we can do the coercion in the ingestor as well as allowing
    it within this pre-checker. Thus, this pre-checker logic will evolve over time as the ingestor
    logic evolves over time.
    """

    # Note: string/string dicts are easier to serialize/deserialize than pa.Schema
    obs_schema: Dict[str, str]
    var_schema: Dict[str, str]
    raw_var_schema: Optional[Dict[str, str]]

    # TODO include 'raw' in X_dtypes or no? Different for AnnData and for SOMA. When in doubt,
    # lean SOMA.
    X_dtypes: Dict[str, str]
    raw_X_dtype: Optional[str]

    obsm_dtypes: Dict[str, str]
    varm_dtypes: Dict[str, str]

    @classmethod
    def from_anndata(
        cls,
        adata: ad.AnnData,
        *,
        default_obs_field_name: str = "obs_id",
        default_var_field_name: str = "var_id",
        default_X_layer_name: str = "data",
    ) -> Self:
        """
        Constructs a pre-check signature from AnnData/H5AD input, which can be compared
        against another signature from AnnData/H5AD or SOMA experiment.

        AnnData inputs have a column offered as the index. This can be: named explicitly "obs_id",
        "var_id", etc.; unnamed: adata.obs.index.name is None; named "index".

        In the latter two cases the ingestor allows a rename to the user's choice such as "obs_id"
        and "var_id". Here in the appender pre-check logic, we allow the same.
        """

        obs_schema = _string_dict_from_pandas_dataframe(
            adata.obs, default_obs_field_name
        )
        var_schema = _string_dict_from_pandas_dataframe(
            adata.var, default_var_field_name
        )

        X_dtypes = {}
        X_dtypes[default_X_layer_name] = str(
            tiledbsoma._arrow_types.arrow_type_from_tiledb_dtype(adata.X.dtype)
        )
        for X_layer_name, X_layer in adata.layers.items():
            X_dtypes[X_layer_name] = str(
                tiledbsoma._arrow_types.arrow_type_from_tiledb_dtype(X_layer.dtype)
            )

        raw_X_dtype = None
        raw_var_schema = None
        if adata.raw is not None:
            raw_X_dtype = str(
                tiledbsoma._arrow_types.arrow_type_from_tiledb_dtype(adata.raw.X.dtype)
            )
            raw_var_schema = _string_dict_from_pandas_dataframe(
                adata.raw.var, default_var_field_name
            )

        obsm_dtypes = {
            k: str(tiledbsoma._arrow_types.arrow_type_from_tiledb_dtype(v.dtype))
            for k, v in adata.obsm.items()
        }
        varm_dtypes = {
            k: str(tiledbsoma._arrow_types.arrow_type_from_tiledb_dtype(v.dtype))
            for k, v in adata.varm.items()
        }

        return cls(
            obs_schema=obs_schema,
            var_schema=var_schema,
            X_dtypes=X_dtypes,
            raw_X_dtype=raw_X_dtype,
            raw_var_schema=raw_var_schema,
            obsm_dtypes=obsm_dtypes,
            varm_dtypes=varm_dtypes,
        )

    @classmethod
    def from_h5ad(
        cls,
        h5ad_file_name: str,
        *,
        default_obs_field_name: str = "obs_id",
        default_var_field_name: str = "var_id",
        default_X_layer_name: str = "data",
    ) -> Self:
        """
        See ``from_anndata``.
        """
        with read_h5ad(h5ad_file_name, mode="r") as adata:
            return cls.from_anndata(
                adata,
                default_X_layer_name=default_X_layer_name,
                default_obs_field_name=default_obs_field_name,
                default_var_field_name=default_var_field_name,
            )

    @classmethod
    def from_soma_experiment(
        cls,
        uri: str,
        measurement_name: str = "RNA",
        context: Optional[SOMATileDBContext] = None,
    ) -> Self:
        """
        Constructs a pre-check signature from a SOMA experiment, which can be compared against
        another signature from AnnData/H5AD or SOMA experiment.
        """

        with tiledbsoma.Experiment.open(uri, context=context) as exp:
            obs_schema = _string_dict_from_arrow_schema(exp.obs.schema)

            var_schema = _string_dict_from_arrow_schema(
                exp.ms[measurement_name].var.schema
            )

            X_dtypes = {}
            for X_layer_name in exp.ms[measurement_name].X.keys():
                X = exp.ms[measurement_name].X[X_layer_name]
                X_dtypes[X_layer_name] = str(X.schema.field("soma_data").type)

            raw_X_dtype = None
            raw_var_schema = None
            if "raw" in exp.ms:
                raw_var_schema = _string_dict_from_arrow_schema(
                    exp.ms["raw"].var.schema
                )

                X = exp.ms["raw"].X[X_layer_name]
                raw_X_dtype = str(X.schema.field("soma_data").type)

            obsm_dtypes: Dict[str, str] = {}
            if "obsm" in exp.ms[measurement_name]:
                for obsm_layer_name in exp.ms[measurement_name].obsm.keys():
                    obsm = exp.ms[measurement_name].obsm[obsm_layer_name]
                    obsm_dtypes[obsm_layer_name] = str(
                        obsm.schema.field("soma_data").type
                    )

            varm_dtypes: Dict[str, str] = {}
            if "varm" in exp.ms[measurement_name]:
                for varm_layer_name in exp.ms[measurement_name].varm.keys():
                    varm = exp.ms[measurement_name].varm[varm_layer_name]
                    varm_dtypes[varm_layer_name] = str(
                        varm.schema.field("soma_data").type
                    )

            return cls(
                obs_schema=obs_schema,
                var_schema=var_schema,
                X_dtypes=X_dtypes,
                raw_X_dtype=raw_X_dtype,
                raw_var_schema=raw_var_schema,
                obsm_dtypes=obsm_dtypes,
                varm_dtypes=varm_dtypes,
            )

    @classmethod
    def check_compatible(cls, signatures: Dict[str, Self]) -> None:
        """
        Determines if any number of signatures from SOMA experiment or AnnData/H5AD will be safe
        from schema-incompatibility at ingestion time. On failure, a ``ValueError`` is raised
        with user-suitable cause details.
        """
        if len(signatures) < 2:
            return
        items = list(signatures.items())
        namea, siga = items[0]
        for nameb, sigb in items[1:]:
            if not siga._compatible_with(sigb):
                raise ValueError(
                    f"Incompatible signatures {namea!r}, {nameb!r}:\n{siga.to_json()}\n{sigb.to_json()}"
                )

    def _compatible_with(self, other: Self) -> bool:
        """
        Pairwise helper method for ``compatible``. Reasons for incompatibility are currently advised
        to be handled a level up by simply showing the user the failed signature pair.
        """

        # Implementation note: _at present_ this could be implemented as `self == other`.  But
        # "coming soon" we'll allow people the ability to do things like coercing one input's
        # float64 to another's float32 and the evolution will be toward more iffing, not less.  As
        # well, in that case, we'll need to also fine-grain the error-reporting to clearly spell out
        # what we coudn't handle -- rather than just showing the user the two signatures' .to_json()
        # output as we do today.

        if self.obs_schema != other.obs_schema:
            return False

        if self.var_schema != other.var_schema:
            return False

        if self.X_dtypes != other.X_dtypes:
            return False

        if self.raw_X_dtype != other.raw_X_dtype:
            return False
        if self.raw_var_schema != other.raw_var_schema:
            return False

        if self.obsm_dtypes != other.obsm_dtypes:
            return False
        if self.varm_dtypes != other.varm_dtypes:
            return False

        return True

    def to_json(self) -> str:
        """Presents a signature as JSON which is suitable for distributed logging."""
        return json.dumps(self, default=attrs.asdict, sort_keys=True, indent=4)

    @classmethod
    def from_json(cls, s: str) -> Self:
        dikt = json.loads(s)
        return cls(**dikt)
