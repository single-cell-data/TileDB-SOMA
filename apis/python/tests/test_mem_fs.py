import numpy as np
import pandas as pd
import pytest
import tiledb
from anndata import AnnData

import tiledbsoma.io as io
from tiledbsoma import SOMA


# TODO: re-enable when #46 is resolved
@pytest.mark.skip(reason="Fails: filed as issue #46")
def test_from_anndata_memfs():
    n_obs = 1000
    n_var = 1
    obs = pd.DataFrame(data={"A": np.arange(n_obs)}, index=np.arange(n_obs).astype(str))
    var = pd.DataFrame(data={"A": np.arange(n_var)}, index=np.arange(n_var).astype(str))
    X = np.ones((n_obs, n_var), dtype=np.float32)
    adata = AnnData(X=X, obs=obs, var=var, dtype=X.dtype)

    path = "mem://foo"
    soma = SOMA(path)
    io.from_anndata(soma, adata)

    with tiledb.Group(soma.uri) as G:
        assert G is not None
        members = [mbr.uri for mbr in G]
        assert f"{path}/X" in members
        assert f"{path}/obs" in members
        assert f"{path}/var" in members
