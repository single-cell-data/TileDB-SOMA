from anndata import AnnData
import tiledb
from tiledbsc import SOMA
import pandas as pd
import numpy as np
from scipy import sparse

from collections import OrderedDict
import os
import pytest

"""
Verify that the AnnData.uns persists correctly. Currently focused on a simple
subset used by cellxgene schema.

TODO: round-trip via to_anndata() when that is implemented.

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
        "float": 3.25,
        "string": "a string",
        "list_of_int": list(i*10 for i in range(10)),
        "list_of_float": list(i*1.25 for i in range(10)),
        "list_of_string": list(str(i*100) for i in range(10)),
        "simple_dict": {"A": 0, "B": "one"},
        "numpy_ndarray_1d_int": np.asarray([1,2,3]),
        "numpy_ndarray_2d_float": np.asarray([[1.,2.,3.],[4.,5.,6.]]),
        "numpy_ndarray_1d_string": np.asarray(['a','b','c']),
        "pandas_dataframe": pd.DataFrame(
            index=np.arange(10).astype(str), data={"A": np.arange(10, dtype=np.int32)}
        ),
    }

    adata = AnnData(X=X, obs=obs, var=var, uns=uns)

    SOMA(tmp_path.as_posix()).from_anndata(adata)

    unspath = tmp_path / "uns"
    assert os.path.exists(unspath)
    for key in uns.keys():
        assert os.path.exists(os.path.join(unspath, key))

    with tiledb.open(os.path.join(unspath, "int")) as A:
        df = A[:]
        assert isinstance(df, np.ndarray)
        assert(df.shape == (1,))
        # Python int goes to int32 by default on Windows
        assert(df.dtype == np.int64 or df.dtype == np.int32)
        assert df[0] == 1

    with tiledb.open(os.path.join(unspath, "float")) as A:
        df = A[:]
        assert isinstance(df, np.ndarray)
        assert(df.shape == (1,))
        assert(df.dtype == np.float64)
        assert df[0] == 3.25

    with tiledb.open(os.path.join(unspath, "string")) as A:
        df = A[:]
        assert isinstance(df, np.ndarray)
        assert(df.shape == (1,))
        assert(df.dtype == np.dtype('O'))
        assert df[0] == "a string"

    with tiledb.open(os.path.join(unspath, "list_of_int")) as A:
        df = A[:]
        assert isinstance(df, np.ndarray)
        assert(df.shape == (10,))
        assert(df.dtype == np.int64 or df.dtype == np.int32)
        assert df[9] == 90

    with tiledb.open(os.path.join(unspath, "list_of_float")) as A:
        df = A[:]
        assert isinstance(df, np.ndarray)
        assert(df.shape == (10,))
        assert(df.dtype == np.float64)
        assert df[9] == 11.25

    with tiledb.open(os.path.join(unspath, "list_of_string")) as A:
        df = A[:]
        assert isinstance(df, np.ndarray)
        assert(df.shape == (10,))
        assert(df.dtype == np.dtype('O'))
        assert df[9] == "900"

    with tiledb.open(os.path.join(unspath, "simple_dict", "A")) as A:
        df = A[:]
        assert isinstance(df, np.ndarray)
        assert(df.shape == (1,))
        assert(df.dtype == np.int64 or df.dtype == np.int32)
        assert df[0] == 0

    with tiledb.open(os.path.join(unspath, "simple_dict", "B")) as A:
        df = A[:]
        assert isinstance(df, np.ndarray)
        assert(df.shape == (1,))
        assert(df.dtype == np.dtype('O'))
        assert df[0] == "one"

    with tiledb.open(os.path.join(unspath, "numpy_ndarray_1d_int")) as A:
        df = A[:]
        assert isinstance(df, np.ndarray)
        assert(df.shape == (3,))
        assert(df.dtype == np.int64 or df.dtype == np.int32)
        assert df[2] == 3

    with tiledb.open(os.path.join(unspath, "numpy_ndarray_2d_float")) as A:
        df = A[:]
        assert isinstance(df, np.ndarray)
        assert(df.shape == (2,3))
        assert(df.dtype == np.float64)
        assert df[1][2] == 6.0

    with tiledb.open(os.path.join(unspath, "numpy_ndarray_1d_string")) as A:
        df = A[:]
        assert isinstance(df, np.ndarray)
        assert(df.shape == (3,))
        assert(df.dtype == np.dtype('O'))
        assert df[2] == 'c'

    with tiledb.open(os.path.join(unspath, "pandas_dataframe")) as A:
        df = A[:]
        assert isinstance(df, OrderedDict)
        dfa = df['A']
        assert isinstance(dfa, np.ndarray)
        assert(dfa.shape == (10,))
        assert(dfa.dtype == np.int32)
        assert dfa[9] == 9
