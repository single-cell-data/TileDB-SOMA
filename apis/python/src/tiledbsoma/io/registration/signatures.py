from dataclasses import dataclass
from typing import Dict, Optional, Sequence, TypeVar

import anndata as ad
import pandas as pd
import numpy as np
import pandas._typing as pdt
from typing_extensions import Self

import tiledbsoma
import tiledbsoma.logging

_DT = TypeVar("_DT", bound=pdt.Dtype)

# ================================================================
# import tiledbsoma.io.registration as reg
# reg.Signature.fromH5AD('/var/s/a/pbmc-small.h5ad')
# ================================================================

@dataclass
class Signature:
    """TODO: docstring"""

    # XXX use arrow -- ?
    obs_index_field_name: str
    obs_schema: Dict[str, _DT]

    var_index_field_name: str
    var_schema: Dict[str, _DT]
    raw_var_schema: Optional[Dict[str, _DT]]

    X_dtypes: Dict[str, _DT]
    raw_X_dtype: Optional[_DT]

    obsm_dtypes: Dict[str, _DT]
    varm_dtypes: Dict[str, _DT]

    @classmethod
    def fromAnnData(cls, adata: ad.AnnData, *, default_X_layer_name: str = "data") -> Self:

        obs_index_field_name = adata.obs.index.name
        obs_schema = adata.obs.dtypes

        var_index_field_name = adata.var.index.name
        var_schema = adata.var.dtypes

        X_dtypes = {}
        X_dtypes[default_X_layer_name] = adata.X.dtype
        for X_layer_name, X_layer in adata.layers.items():
            X_dtypes[X_layer_name] = X_layer.dtype
            
        raw_X_dtype = None if adata.raw is None else adata.raw.X.dtype
        raw_var_schema = None if adata.raw is None else adata.raw.var.dtypes

        obsm_dtypes = { k:v.dtype for k,v in adata.obsm.items() }
        varm_dtypes = { k:v.dtype for k,v in adata.varm.items() }

        return cls(
            obs_index_field_name=obs_index_field_name,
            obs_schema=obs_schema,
            var_index_field_name=var_index_field_name,
            var_schema=var_schema,
            X_dtypes=X_dtypes,
            raw_X_dtype=raw_X_dtype,
            raw_var_schema=raw_var_schema,
            obsm_dtypes=obsm_dtypes,
            varm_dtypes=varm_dtypes,
        )

    @classmethod
    def fromH5AD(cls, h5ad_file_name: str, *, default_X_layer_name: str = "data") -> Self:
        adata = ad.read_h5ad(h5ad_file_name, "r")
        return cls.fromAnnData(adata, default_X_layer_name=default_X_layer_name)

    @classmethod
    def fromSOMAExperiment(cls, uri: str, measurement_name: str = "RNA") -> Self:
        with tiledbsoma.Experiment.open(uri) as exp:

            # The coord 0 is arbitrary. We don't need to read _any_ data -- we just need the schema.
            # The SOMA API has a no-data-read Arrow-schema accessors, but we want the Pandas/NumPy
            # schema to check against AnnData inputs.
            obs_df = exp.obs.read(coords=[0]).concat().to_pandas()
            obs_index_field_name = obs_df.index.name
            obs_schema = { k:v.dtype for k,v in obs_df.items() }

            var_df = exp.ms[measurement_name].var.read(coords=[0]).concat().to_pandas()
            var_index_field_name = var_df.index.name
            var_schema = { k:v.dtype for k,v in var_df.items() }

            X_dtypes = {}
            for X_layer_name in exp.ms[measurement_name].X.keys():
                X = exp.ms[measurement_name].X[X_layer_name].read(coords=[0,0]).to_numpy()
                X_dtypes[X_layer_name] = X.dtype

            raw_X_dtype = None
            raw_var_schema = None
            # TODO

            # TODO
            obsm_dtypes = {}
            varm_dtypes = {}

            return cls(
                obs_index_field_name=obs_index_field_name,
                obs_schema=obs_schema,
                var_index_field_name=var_index_field_name,
                var_schema=var_schema,
                X_dtypes=X_dtypes,
                raw_X_dtype=raw_X_dtype,
                raw_var_schema=raw_var_schema,
                obsm_dtypes=obsm_dtypes,
                varm_dtypes=varm_dtypes,
            )

    @classmethod
    def compatible(cls, signatures: Sequence[Self]) -> bool:
        if len(signatures) < 2:
            return True
        for other in signatures[1:]:
            if not signatures[0].compatibleWith(other):
                return False
        return True

    # TODO: return why-not as loggable message, and/or log it directly here
    def compatibleWith(self, other: Self) -> bool:
        # TODO: if None, handleable
        if self.obs_index_field_name != other.obs_index_field_name:
            return False

        if any(self.obs_schema != other.obs_schema):
            return False

        if self.var_index_field_name != other.var_index_field_name:
            return False

        if any(self.var_schema != other.var_schema):
            return False

        # TODO: more

        if self.X_dtypes != other.X_dtypes:
            return False

        if self.raw_X_dtype != other.raw_X_dtype:
            return False
        if any(self.raw_var_schema != other.raw_var_schema):
            return False

        if self.obsm_dtypes != other.obsm_dtypes:
            return False
        if self.varm_dtypes != other.varm_dtypes:
            return False

        return True
