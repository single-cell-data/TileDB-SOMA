from typing import Dict, Tuple

import anndata
import pyarrow as pa
import scipy.sparse as sparse

from .types import ExperimentQueryReadArrowResult


def arrow_to_scipy_csr(X: pa.Table, shape: Tuple[int, int]) -> sparse.csr_matrix:
    return sparse.csr_matrix(
        (X["soma_data"].to_numpy(), (X["_dim_0"].to_numpy(), X["_dim_1"].to_numpy())),
        shape=shape,
    )


def make_anndata(query_result: ExperimentQueryReadArrowResult) -> anndata.AnnData:

    obs = query_result["obs"]
    obs = obs.to_pandas()
    obs.index = obs.index.map(str)

    var = query_result["var"]
    var = var.to_pandas()
    var.index = var.index.map(str)

    shape = (len(obs), len(var))

    X = query_result.get("X", None)
    if X is not None:
        X = arrow_to_scipy_csr(X, shape)

    X_layers = query_result.get("X_layers", {})
    layers: Dict[str, sparse.csr_matrix] = {}
    for X_layer_name, X_layer_table in X_layers.items():
        layers[X_layer_name] = arrow_to_scipy_csr(X_layer_table, shape)

    return anndata.AnnData(
        X=X, obs=obs, var=var, layers=(layers if len(layers) else None)
    )
