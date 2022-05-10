from anndata import AnnData
import tiledb
from tiledbsc import SOMA
import pandas as pd
import numpy as np
from scipy import sparse
import os

import pytest

"""
Verify that the AnnData.uns persists correctly. Currently focused on a simple
subset used by cellxgene schema.

TODO: round-trip via to_anndata() when that is implemented.
TODO: enable test when `uns` supported (issue #54)

"""

def test_from_anndata_uns(tmp_path):
    """
    Test very simple `uns` is persisted.
    """

    X = np.ones((10, 10), dtype=np.float32)
    obs = pd.DataFrame(
          index=np.arange(10).astype(str), data={"A": np.arange(10, dtype=np.int32)}
    )
    var = pd.DataFrame(
          index=np.arange(10).astype(str), data={"A": np.arange(10, dtype=np.int32)}
    )

    uns = {
        "int": 1,
        "float": 2.3,
        "string": "a string",
        "list_of_string": list(str(i) for i in range(10)),
        "list_of_int": list(i*10 for i in range(10)),
        "list_of_float": list(i*1.25 for i in range(10)),
        "simple_dict": {"A": 0, "B": "one"},
        "numpy_ndarray_1d_int": np.asarray([1,2,3]),
        "numpy_ndarray_2d_float": np.asarray([[1.,2.,3.],[4.,5.,6.]]),
        "numpy_ndarray_1d_string": np.asarray(['a','b','c']),
        "pandas_dataframe": pd.DataFrame(
            index=np.arange(10).astype(str), data={"A": np.arange(10, dtype=np.int32)}
        ),
    }

# int
# float
# list_of_float
# list_of_int
# list_of_string
# numpy_ndarray_1d_int
# numpy_ndarray_1d_string
# numpy_ndarray_2d_float
# pandas_dataframe
# simple_dict
# string

    adata = AnnData(X=X, obs=obs, var=var, uns=uns)
    SOMA(tmp_path.as_posix()).from_anndata(adata)

    # TODO: persistence for UNS is not yet defined. Update when it is.
    assert (tmp_path / "uns").exists()
    for key in uns.keys():
        assert (tmp_path / "uns" / key).exists()
        # TODO: add check that contents were saved correctly
