import json
from typing import Dict, Optional

import anndata as ad
import attrs
import pandas as pd
import pyarrow as pa
from typing_extensions import Self

import tiledbsoma
import tiledbsoma.logging
from tiledbsoma._arrow_types import df_to_arrow
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
    Here we provide compatiblity with the ingestor.

    SOMA experiments are indexed by int64 soma_joinid and this is SOMA-only.

    AnnData inputs have a column offered as the index. This can be: named explicitly "obs_id",
    "var_id", etc.; unnamed: adata.obs.index.name is None; named "index".

    In the latter two cases the ingestor allows a rename to the user's choice
    such as "obs_id" and "var_id". Here in the appender pre-check logic, we
    allow the same.
    """

    df = df.head(1)  # since reset_index can be expensive on full data
    if df.index.name is None or df.index.name == "index":
        df.reset_index(inplace=True)
        if default_index_name in df:
            if "index" in df:
                # Avoid the warning:
                # "A value is trying to be set on a copy of a slice from a DataFrame"
                # which would occur if we did:
                # df.drop(columns=["index"], inplace=True)
                df = df.drop(columns=["index"])
        else:
            df.rename(columns={"index": default_index_name}, inplace=True)
    else:
        df.reset_index(inplace=True)

    arrow_table = df_to_arrow(df)
    arrow_schema = arrow_table.schema.remove_metadata()
    return _string_dict_from_arrow_schema(arrow_schema)


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
        adata = ad.read_h5ad(h5ad_file_name, "r")
        return cls.from_anndata(adata, default_X_layer_name=default_X_layer_name)

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
            obsm_dtypes = {}
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
