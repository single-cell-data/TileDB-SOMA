import json
from dataclasses import dataclass
from typing import Dict, Optional, Tuple

import anndata as ad
import pandas as pd
import pyarrow as pa
from typing_extensions import Self

import tiledbsoma
import tiledbsoma.logging


# TODO TEMP COPY FROM VIO
def _df_to_arrow(df: pd.DataFrame) -> pa.Table:
    """
    Categoricals are not yet well supported, so we must flatten.
    Also replace Numpy/Pandas-style nulls with Arrow-style nulls.
    """
    null_fields = set()
    for k in df:
        if df[k].dtype == "category":
            df[k] = df[k].astype(df[k].cat.categories.dtype)
        if df[k].isnull().any():
            if df[k].isnull().all():
                df[k] = pa.nulls(df.shape[0], pa.infer_type(df[k]))
            else:
                df[k].where(
                    df[k].notnull(),
                    pd.Series(pa.nulls(df[k].isnull().sum(), pa.infer_type(df[k]))),
                    inplace=True,
                )
            null_fields.add(k)
    arrow_table = pa.Table.from_pandas(df)
    if null_fields:
        md = arrow_table.schema.metadata
        md.update(dict.fromkeys(null_fields, "nullable"))
        arrow_table = arrow_table.replace_schema_metadata(md)
    return arrow_table


def _string_dict_from_arrow_schema(schema: pa.Schema) -> Dict[str, str]:
    """
    Converts an Arrow schema to a string/string dict, which is easier on the eyes,
    easier to convert from/to JSON, and easier to do del-key on.
    """

    def stringify_type(t: pa.DataType) -> str:
        retval = str(t)
        if retval == "large_string":
            return "string"  # TODO comment
        if retval == "large_binary":
            return "binary"  # TODO comment
        return retval

    retval = {name: stringify_type(schema.field(name).type) for name in schema.names}
    if "soma_joinid" in retval:
        del retval["soma_joinid"]  # TODO comment
    return retval


def _string_dict_from_pandas_dataframe(
    df: pd.DataFrame,
    default_index_name: str,
) -> Dict[str, str]:
    df = df.head(1)  # since reset_index can be expensive on full data
    # TODO: comment relative to vio
    if df.index.name is None or df.index.name == "index":
        df.reset_index(inplace=True)
        df.rename(columns={"index": default_index_name}, inplace=True)
    else:
        df.reset_index(inplace=True)

    arrow_table = _df_to_arrow(df)
    arrow_schema = arrow_table.schema.remove_metadata()
    return _string_dict_from_arrow_schema(arrow_schema)


@dataclass
class Signature:
    """TODO: docstring"""

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
    def fromAnnData(
        cls,
        adata: ad.AnnData,
        *,
        default_obs_field_name: str = "obs_id",  # TODO: describe
        default_var_field_name: str = "var_id",  # TODO: describe
        default_X_layer_name: str = "data",
    ) -> Self:

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
    def fromH5AD(
        cls,
        h5ad_file_name: str,
        *,
        default_obs_field_name: str = "obs_id",
        default_var_field_name: str = "var_id",
        default_X_layer_name: str = "data",
    ) -> Self:
        adata = ad.read_h5ad(h5ad_file_name, "r")
        return cls.fromAnnData(adata, default_X_layer_name=default_X_layer_name)

    @classmethod
    def fromSOMAExperiment(cls, uri: str, measurement_name: str = "RNA") -> Self:
        with tiledbsoma.Experiment.open(uri) as exp:

            obs_schema = _string_dict_from_arrow_schema(exp.obs.schema)

            var_schema = _string_dict_from_arrow_schema(
                exp.ms[measurement_name].var.schema
            )

            X_dtypes = {}
            for X_layer_name in exp.ms[measurement_name].X.keys():
                X = exp.ms[measurement_name].X[X_layer_name]
                # TODO: replace [2] with right lookup for the soma_data attr
                X_dtypes[X_layer_name] = str(X.schema[2].type)

            raw_X_dtype = None
            raw_var_schema = None
            if "raw" in exp.ms:
                raw_var_schema = _string_dict_from_arrow_schema(
                    exp.ms["raw"].var.schema
                )

                X = exp.ms["raw"].X[X_layer_name]
                # TODO: replace [2] with right lookup for the soma_data attr
                raw_X_dtype = str(X.schema[2].type)

            obsm_dtypes: Dict[str, str] = {}
            obsm_dtypes = {}
            if "obsm" in exp.ms[measurement_name]:
                for obsm_layer_name in exp.ms[measurement_name].obsm.keys():
                    obsm = exp.ms[measurement_name].obsm[obsm_layer_name]
                    # TODO: replace [2] with right lookup for the soma_data attr
                    obsm_dtypes[obsm_layer_name] = str(obsm.schema[2].type)

            varm_dtypes: Dict[str, str] = {}
            if "varm" in exp.ms[measurement_name]:
                for varm_layer_name in exp.ms[measurement_name].varm.keys():
                    varm = exp.ms[measurement_name].varm[varm_layer_name]
                    # TODO: replace [2] with right lookup for the soma_data attr
                    varm_dtypes[varm_layer_name] = str(varm.schema[2].type)

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
    def compatible(cls, signatures: Dict[str, Self]) -> Tuple[bool, Optional[str]]:
        if len(signatures) < 2:
            return (True, None)
        names = list(signatures.keys())
        first_name = names[0]
        for other_name in names[1:]:
            siga = signatures[first_name]
            sigb = signatures[other_name]
            if not siga.compatibleWith(sigb):
                msg = f"Incompatible signatures {first_name!r}, {other_name!r}:\n{siga.toJSON()}\n{sigb.toJSON()}"
                return (False, msg)
        return (True, None)

    # TODO: return why-not as loggable message(s), and/or log it/them directly here
    # Or: maybe it clear folks can use .toJSON() and inspect ...
    def compatibleWith(self, other: Self) -> bool:

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

    def toJSON(self) -> str:
        """TODO: docstring"""
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, indent=4)

    @classmethod
    def fromJSON(cls, s: str) -> Self:
        dikt = json.loads(s)
        return cls(
            dikt["obs_schema"],
            dikt["var_schema"],
            dikt["raw_var_schema"],
            dikt["X_dtypes"],
            dikt["raw_X_dtype"],
            dikt["obsm_dtypes"],
            dikt["varm_dtypes"],
        )
