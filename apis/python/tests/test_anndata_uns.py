import os
from collections import OrderedDict

import numpy as np
import pandas as pd
import tiledb
from anndata import AnnData

import tiledbsoma.io as io
from tiledbsoma import SOMA

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
        "scalar_int": 1,
        "scalar_float": 3.25,
        "scalar_string": "a string",
        "list_of_int": list(i * 10 for i in range(10)),
        "list_of_float": list(i * 1.25 for i in range(10)),
        "list_of_string": list(str(i * 100) for i in range(10)),
        "simple_dict": {"A": 0, "B": "one"},
        "numpy_ndarray_1d_int": np.asarray([1, 2, 3]),
        "numpy_ndarray_2d_float": np.asarray([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]),
        "numpy_ndarray_1d_string": np.asarray(["a", "b", "c"]),
        "pandas_dataframe": pd.DataFrame(
            index=np.arange(10).astype(str), data={"A": np.arange(10, dtype=np.int32)}
        ),
    }

    adata = AnnData(X=X, obs=obs, var=var, uns=uns)

    io.from_anndata(SOMA(tmp_path.as_posix()), adata)

    # Example of what we're verifying:

    # >>> soma.uns
    #
    # uns:
    # scalar_float: 3.25
    # scalar_string: a string
    # uns/scalar_int/
    # [1]
    # uns/list_of_float/
    # [ 0.    1.25  2.5   3.75  5.    6.25  7.5   8.75 10.   11.25]
    # uns/list_of_int/
    # [ 0 10 20 30 40 50 60 70 80 90]
    # uns/numpy_ndarray_1d_int/
    # [1 2 3]
    # uns/numpy_ndarray_2d_float/
    # [[1. 2. 3.]
    #  [4. 5. 6.]]
    #
    # uns/simple_dict:
    # B: one
    # uns/simple_dict/A/
    # [0]
    # uns/numpy_ndarray_1d_string/
    # ['a' 'b' 'c']
    # uns/pandas_dataframe/
    # OrderedDict([('A', array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], dtype=int32)), ('__tiledb_rows', array([b'0', b'1', b'2', b'3', b'4', b'5', b'6', b'7', b'8', b'9'],
    #       dtype=object))])
    # uns/list_of_string/
    # ['0' '100' '200' '300' '400' '500' '600' '700' '800' '900']

    unspath = tmp_path / "uns"
    assert os.path.exists(unspath)
    for key in uns.keys():
        if not key.startswith("scalar_"):
            # Scalars are written as group metadata
            assert (unspath / key).exists()

        with tiledb.group.Group(unspath.as_posix()) as G:
            assert G.meta["scalar_int"] == 1
            assert G.meta["scalar_float"] == 3.25
            assert G.meta["scalar_string"] == "a string"

    with tiledb.open((unspath / "list_of_int").as_posix()) as A:
        df = A[:]
        assert isinstance(df, np.ndarray)
        assert df.shape == (10,)
        assert df.dtype == np.int64 or df.dtype == np.int32
        assert df[9] == 90

    with tiledb.open((unspath / "list_of_float").as_posix()) as A:
        df = A[:]
        assert isinstance(df, np.ndarray)
        assert df.shape == (10,)
        assert df.dtype == np.float64
        assert df[9] == 11.25

    with tiledb.open((unspath / "list_of_string").as_posix()) as A:
        df = A[:]
        assert isinstance(df, np.ndarray)
        assert df.shape == (10,)
        assert df.dtype == np.dtype("O")
        assert df[9] == "900"

    with tiledb.group.Group((unspath / "simple_dict").as_posix()) as G:
        # Scalars are written as group metadata
        assert G.meta["A"] == 0
        assert G.meta["B"] == "one"

    with tiledb.open((unspath / "numpy_ndarray_1d_int").as_posix()) as A:
        df = A[:]
        assert isinstance(df, np.ndarray)
        assert df.shape == (3,)
        assert df.dtype == np.int64 or df.dtype == np.int32
        assert df[2] == 3

    with tiledb.open((unspath / "numpy_ndarray_2d_float").as_posix()) as A:
        df = A[:]
        assert isinstance(df, np.ndarray)
        assert df.shape == (2, 3)
        assert df.dtype == np.float64
        assert df[1][2] == 6.0

    with tiledb.open((unspath / "numpy_ndarray_1d_string").as_posix()) as A:
        df = A[:]
        assert isinstance(df, np.ndarray)
        assert df.shape == (3,)
        assert df.dtype == np.dtype("O")
        assert df[2] == "c"

    with tiledb.open((unspath / "pandas_dataframe").as_posix()) as A:
        df = A[:]
        assert isinstance(df, OrderedDict)
        dfa = df["A"]
        assert isinstance(dfa, np.ndarray)
        assert dfa.shape == (10,)
        assert dfa.dtype == np.int32
        assert dfa[9] == 9
