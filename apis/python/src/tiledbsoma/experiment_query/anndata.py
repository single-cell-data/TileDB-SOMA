from typing import Dict, Tuple

import anndata
import pyarrow as pa
import scipy.sparse as sparse

from .eq_types import ExperimentAxisQueryReadArrowResult


def _arrow_to_scipy_csr(X: pa.Table, shape: Tuple[int, int]) -> sparse.csr_matrix:
    """
    Private utility which converts a table repesentation of X to a CSR matrix.

    IMPORTANT: by convention, assumes that the data is positionally indexed (hence
    the use of _dim_{n} rather than soma_dim{n}).

    See query.py::_rewrite_X_for_positional_indexing for more info.
    """
    assert "_dim_0" in X.column_names, "X must be positionally indexed"
    assert "_dim_1" in X.column_names, "X must be positionally indexed"

    return sparse.csr_matrix(
        (X["soma_data"].to_numpy(), (X["_dim_0"].to_numpy(), X["_dim_1"].to_numpy())),
        shape=shape,
    )


def _make_anndata(query_result: ExperimentAxisQueryReadArrowResult) -> anndata.AnnData:
    """Private utility to create an AnnData from an ExperimentAxisQueryReadArrowResult"""

    obs = query_result["obs"]
    obs = obs.to_pandas()
    obs.index = obs.index.map(str)

    var = query_result["var"]
    var = var.to_pandas()
    var.index = var.index.map(str)

    shape = (len(obs), len(var))

    X = query_result.get("X", None)
    if X is not None:
        X = _arrow_to_scipy_csr(X, shape)

    X_layers = query_result.get("X_layers", {})
    layers: Dict[str, sparse.csr_matrix] = {}
    for X_layer_name, X_layer_table in X_layers.items():
        layers[X_layer_name] = _arrow_to_scipy_csr(X_layer_table, shape)

    return anndata.AnnData(
        X=X, obs=obs, var=var, layers=(layers if len(layers) else None)
    )
