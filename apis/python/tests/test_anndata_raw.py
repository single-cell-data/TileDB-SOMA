"""
Tests which verify from_anndata behavior with `.raw`.

https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.raw.html#anndata.AnnData.raw
https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html

Conventionally, `.raw` is a subset of var and X. In practice, this does not appear to be
enforced by AnnData/ScanPy.
"""

import numpy as np
import pandas as pd
import pytest
import tiledb
from anndata import AnnData

import tiledbsoma.io as io
from tiledbsoma import SOMA


def test_from_anndata_raw_X(tmp_path, adata):
    """
    Verify that anndata.raw.X is correctly saved.
    """
    io.from_anndata(SOMA(tmp_path.as_posix()), adata)

    assert all(
        (tmp_path / sub_array_path).exists()
        for sub_array_path in ["X/data", "raw/X/data"]
    )

    with tiledb.open((tmp_path / "raw" / "X" / "data").as_posix()) as raw_X:
        orig_raw_X_df = pd.DataFrame(
            adata.raw.X.flatten(),
            index=pd.MultiIndex.from_product(
                [adata.obs.index, adata.raw.var.index], names=["obs_id", "var_id"]
            ),
            columns=["value"],
        ).reset_index()
        assert raw_X.schema.sparse
        assert raw_X.df[:].equals(orig_raw_X_df)

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


def test_from_anndata_raw_var(tmp_path, adata):
    """
    Verify that anndata.raw.var is correctly saved.
    """
    io.from_anndata(SOMA(tmp_path.as_posix()), adata)

    # TODO: no idea where `.raw.var` will be saved. Update when issue #50 is resolved
    assert all(
        (tmp_path / sub_array_path).exists() for sub_array_path in ["var", "raw/var"]
    )

    with tiledb.open((tmp_path / "raw" / "var").as_posix()) as raw_var:
        assert raw_var.schema.sparse
        assert raw_var.df[:].equals(adata.raw.var)


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
