from typing import Tuple

import anndata
import pyarrow as pa
import scipy.sparse as sparse

from .eq_types import ExperimentAxisQueryReadArrowResult


def _arrow_to_scipy_csr(
    arrow_table: pa.Table, shape: Tuple[int, int]
) -> sparse.csr_matrix:
    """
    Private utility which converts a table repesentation of X to a CSR matrix.

    IMPORTANT: by convention, assumes that the data is positionally indexed (hence
    the use of _dim_{n} rather than soma_dim{n}).

    See query.py::_rewrite_X_for_positional_indexing for more info.
    """
    assert "_dim_0" in arrow_table.column_names, "X must be positionally indexed"
    assert "_dim_1" in arrow_table.column_names, "X must be positionally indexed"

    return sparse.csr_matrix(
        (
            arrow_table["soma_data"].to_numpy(),
            (arrow_table["_dim_0"].to_numpy(), arrow_table["_dim_1"].to_numpy()),
        ),
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

    x = query_result.get("X")
    if x is not None:
        x = _arrow_to_scipy_csr(x, shape)

    x_layers = query_result.get("X_layers", {})
    layers = {
        name: _arrow_to_scipy_csr(table, shape) for name, table in x_layers.items()
    }

    return anndata.AnnData(X=x, obs=obs, var=var, layers=(layers or None))
