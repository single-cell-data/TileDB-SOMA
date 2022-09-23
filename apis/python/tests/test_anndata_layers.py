import numpy as np
import pandas as pd
import pytest
import tiledb
from anndata import AnnData

import tiledbsoma.io as io
from tiledbsoma import SOMA


def test_from_anndata_X_layers(tmp_path, adata):
    adata.layers["extra"] = adata.X

    io.from_anndata(SOMA(tmp_path.as_posix()), adata)

    assert (tmp_path / "X" / "data").exists()
    assert (tmp_path / "X" / "extra").exists()

    with tiledb.open((tmp_path / "X" / "data").as_posix()) as X:
        orig_X_df = pd.DataFrame(
            adata.X.flatten(),
            index=pd.MultiIndex.from_product(
                [adata.obs.index, adata.var.index], names=["obs_id", "var_id"]
            ),
            columns=["value"],
        ).reset_index()
        assert X.schema.sparse
        assert X.df[:].equals(orig_X_df)

    with tiledb.open((tmp_path / "X" / "extra").as_posix()) as X:
        orig_X_df = pd.DataFrame(
            adata.X.flatten(),
            index=pd.MultiIndex.from_product(
                [adata.obs.index, adata.var.index], names=["obs_id", "var_id"]
            ),
            columns=["value"],
        ).reset_index()
        assert X.schema.sparse
        assert X.df[:].equals(orig_X_df)


@pytest.fixture
def adata():
    n_obs = 10
    n_var = 8
    obs = pd.DataFrame(
        index=np.arange(n_obs).astype(str), data={"A": np.arange(n_obs, dtype=np.int32)}
    )
    var = pd.DataFrame(
        index=np.arange(n_var).astype(str), data={"A": np.arange(n_var, dtype=np.int32)}
    )
    # Note: no zeros in X
    X_values = np.arange(1, n_obs * n_var + 1, dtype=np.float32)
    X = X_values.reshape((n_obs, n_var))

    adata = AnnData(X=X, obs=obs, var=var, dtype=X.dtype)

    # following the ScanPy/AnnData recipe for filtering genes
    adata.raw = adata
    adata = adata[:, 1 : n_var // 2 : 2]

    return adata
