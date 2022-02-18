#%%
import itertools
import os
import sys
from dataclasses import dataclass
from typing import List
from unittest import result
import warnings

import anndata as ad
import numpy as np
import pandas as pd
import scanpy
import tiledb

#%%
@dataclass
class AnnSchemaBuilder:
    """Create appropriate schema from annotation frame"""

    ann: pd.DataFrame

    def visit_dims(self):
        if self.ann.index is not None and self.ann.index.name is not None:
            name = self.ann.index.name
        else:
            name = ""
        return [tiledb.Dim(name=name, dtype="ascii")]

    def visit_attrs(self):
        result_attrs = []
        for name, col in self.ann.items():
            dtype = tiledb.dataframe_.ColumnInfo.from_values(col)
            result_attrs.append(tiledb.Attr(name=name, dtype=dtype))

        return result_attrs

    def to_schema(self):
        dims = self.visit_dims()
        attrs = self.visit_attrs()

        return tiledb.ArraySchema(domain=tiledb.Domain(*dims), attrs=attrs, sparse=True)


def create_X_schema(dim1: str, dim2: str, *, data_dtype=np.float32):
    # TODO TBD layers need multiple attributes
    return tiledb.ArraySchema(
        domain=tiledb.Domain(
            tiledb.Dim(name=dim1, dtype="ascii"), tiledb.Dim(name=dim2, dtype="ascii")
        ),
        attrs=[tiledb.Attr(name="data", dtype=data_dtype)],
        sparse=True,
    )


class AnnDataExporter:
    """Write AnnData in-memory to TileDB"""

    data: ad.AnnData

    def __init__(self, anndata):
        self.data = anndata

    @classmethod
    def from_10x(cls, path: str):
        warnings.warn("(tmp) writing subset")
        anndata = scanpy.read_10x_h5(path)[:5000]

        warnings.warn("Calling .var_names_make_uniqe on AnnData object!")
        anndata.var_names_make_unique()

        return cls(anndata)

    @classmethod
    def from_h5ad(cls, path: str):
        warnings.warn("(tmp) writing subset")
        anndata = scanpy.read_h5ad(path)[:5000]

        warnings.warn("Calling .var_names_make_uniqe on AnnData object!")
        anndata.var_names_make_unique()

        # TMP TODO convert from categorical to series
        warnings.warn("Converting from categorical to object series!")
        uncat = lambda x: x.astype("O") if isinstance(x.dtype, pd.CategoricalDtype) else x
        var = pd.DataFrame.from_dict({k: uncat(v) for k,v in anndata.var.items()})

        # Make a new object with the modified data
        anndata = ad.AnnData(X=anndata.X, var=var, obs=anndata.obs, raw=anndata.raw)

        return cls(anndata)


    def _create_arrays(self, base: str):
        obs_schema = AnnSchemaBuilder(self.data.obs).to_schema()
        var_schema = AnnSchemaBuilder(self.data.var).to_schema()

        # TODO how does AnnData id the primary keys for the obs/var dfs?
        # TODO create X schema, use 1st dim names from obs and var if available?
        # obs_schema.domain.dim(0).name,
        # var_schema.domain.dim(0).name
        dim1 = "obs"
        dim2 = "var"

        X_schema = create_X_schema(dim1, dim2)
        tiledb.group_create(str(base))

        self._obs_path = os.path.join(base, "obs")
        tiledb.Array.create(str(self._obs_path), obs_schema)

        self._var_path = os.path.join(base, "var")
        tiledb.Array.create(str(self._var_path), var_schema)

        self._data_path = os.path.join(base, "X")
        tiledb.Array.create(str(self._data_path), X_schema)

    def _write_ann(self, path: str, df: pd.DataFrame):
        coords = df.index  # TODO support multi-index?
        data = {name: data for name, data in df.items()}

        with tiledb.open(str(path), "w") as T:
            T[coords] = data

    def _write_data(self):
        coo = self.data.X.tocoo()
        obs_coord = self.data.obs.iloc[coo.row].index
        var_coord = self.data.var.iloc[coo.col].index

        with tiledb.open(str(self._data_path), "w") as T:
            T[obs_coord, var_coord] = {"data": coo.data}

    def write_tiledb(self, base_uri: str):
        self._create_arrays(base_uri)

        self._write_ann(self._obs_path, self.data.obs)
        self._write_ann(self._var_path, self.data.var)
        self._write_data()

if __name__ == "__main__":
    input_path = sys.argv[1]
    output_path = sys.argv[2]

    if input_path.endswith("h5ad"):
        exporter = AnnDataExporter.from_h5ad(input_path)
    else:
        exporter = AnnDataExporter.from_10x(input_path)

    exporter.write_tiledb(output_path)
