from anndata import AnnData
import tiledb
from tiledbsc import SOMA
import pandas as pd
import numpy as np
from scipy import sparse

import pytest

"""
Testing `from_anndata` with the wide diversity of types latent in AnnData.

Status:
* Currently focuses on X, obs and var
* TODO: obsm, varm, uns, et al.
* TODO: re-enable tests that are disabled due to known issues
"""


X_type_sweep = [
    (dtype_name, encoding)
    for dtype_name in [
        "float16",    # TODO: Enable when #39 is fixed
        "float32",
        "float64",
        "int8",
        "int16",
        "int32",
        "int64",
        "uint8",
        "uint16",
        "uint32",
        "uint64",
    ]
    for encoding in ["dense", "csc", "csr"]
]


@pytest.mark.parametrize("X_dtype_name,X_encoding", X_type_sweep)
def test_from_anndata_X_type(tmp_path, X_dtype_name, X_encoding):
    """
    Verify X matrix converts with integrity:
    * X array is created and has correct cardinality
    * X array type `kind` is the same as the source
    * X array type does not lose information/precision
    * X array has expected sparsity
    """
    n_obs = 100
    n_var = 1
    obs = pd.DataFrame(data={"A": np.arange(n_obs, dtype=np.int32)})
    var = pd.DataFrame(data={"A": np.arange(n_var, dtype=np.int32)})

    X_dtype = np.dtype(X_dtype_name)
    if X_encoding == "dense":
        X = np.eye(n_obs, n_var, dtype=X_dtype)
    elif X_encoding == "csc":
        X = sparse.eye(n_obs, n_var, dtype=X_dtype).tocsc()
    elif X_encoding == "csr":
        X = sparse.eye(n_obs, n_var, dtype=X_dtype).tocsr()
    else:
        assert False  # sanity - test misconfiguration

    adata = AnnData(X=X, obs=obs, var=var, dtype=X.dtype)
    print(" =============================================================>==", adata.X.dtype, X_dtype)
    #if X_dtype == np.float16:
        #print(" ????=========================================================>==", adata.X.dtype, X_dtype)
        #assert adata.X.dtype == np.float32
    #else:
        #assert adata.X.dtype == X_dtype  # sanity
    assert adata.X.dtype == X_dtype  # sanity

    SOMA(tmp_path.as_posix()).from_anndata(adata)
    assert all(
        (tmp_path / sub_array_path).exists()
        for sub_array_path in ["obs", "var", "X/data"]
    )

    # check X type & shape.
    with tiledb.open((tmp_path / "X" / "data").as_posix()) as X:
        assert adata.X.dtype.kind == X.schema.attr(0).dtype.kind
        assert adata.X.dtype <= X.schema.attr(0).dtype

        # TODO: at some point, we should be able to posit and validate the
        # sparsity of the resulting X array. At the moment, the check is
        # agnostic of the source sparsity.

        if X.schema.sparse:
            if hasattr(adata.X, "nnz"):
                assert adata.X.nnz == len(X.query(dims=[]).df[:])
            else:
                assert np.count_nonzero(adata.X) == len(X.query(dims=[]).df[:])

        else:
            if hasattr(adata.X, "nnz"):
                assert adata.X.nnz == np.np.count_nonzero(
                    X.query(dims=[]).multi_index[:]["value"]
                )
            else:
                assert adata.X.size == X.query(dims=[]).multi_index[:]["value"]


def test_from_anndata_DataFrame_type(tmp_path):
    """
    Sweep all types that we may see in a DataFrame, and ensure if converts with integrity.
    """
    n = 10
    df_col_type_sweep = [
        ("bool", lambda a: a.astype(bool)),
        ("str", lambda a: a.astype(str)),
        ("bytes", lambda a: a.astype(str).astype(bytes)),
        # ("float16", lambda a: a.astype(np.dtype("float16"))),         TODO: Enable when #39 is fixed
        ("float32", lambda a: a.astype("float32")),
        ("float64", lambda a: a.astype("float64")),
        ("int8", lambda a: a.astype("int8")),
        ("int16", lambda a: a.astype("int16")),
        ("int32", lambda a: a.astype("int32")),
        ("int64", lambda a: a.astype("int64")),
        ("uint8", lambda a: a.astype("uint8")),
        ("uint16", lambda a: a.astype("uint16")),
        ("uint32", lambda a: a.astype("uint32")),
        ("uint64", lambda a: a.astype("uint64")),
        ("object", lambda a: a.astype(str).astype(np.dtype("O"))),
        ("categorical(str)", lambda a: a.astype(str).astype("category")),
        # TODO: The following tests fail due to issue #30 -- re-enable test when resolved.
        # (
        #     "categorical(int32)",
        #     lambda a: a.astype("int32").astype(
        #         pd.CategoricalDtype(categories=a.astype("int32"))
        #     ),
        # ),
        # (
        #     "categorical(uint32)",
        #     lambda a: a.astype("uint32").astype(
        #         pd.CategoricalDtype(categories=a.astype("uint32"))
        #     ),
        # ),
        # (
        #     "categorical(float32)",
        #     lambda a: a.astype("float32").astype(
        #         pd.CategoricalDtype(categories=a.astype("float32"))
        #     ),
        # ),
        # (
        #     "categorical(bool)",
        #     lambda a: a.astype("bool").astype(
        #         pd.CategoricalDtype(categories=a.astype("bool").unique())
        #     ),
        # ),
    ]
    index = np.arange(n).astype(str)  # AnnData requires string indices
    df = pd.DataFrame(
        data={
            f"col_{name}": cast(pd.Series(index=index, data=np.arange(n)))
            for name, cast in df_col_type_sweep
        },
    )
    X = np.ones((n, n), dtype=np.float32)
    adata = AnnData(X=X, obs=df, var=df, dtype=X.dtype)
    SOMA(tmp_path.as_posix()).from_anndata(adata)
    assert all(
        (tmp_path / sub_array_path).exists()
        for sub_array_path in ["obs", "var", "X/data"]
    )

    def cmp_dtype(series, tdb: tiledb.Attr) -> bool:
        """Encapsulate expected conversions moving into TileDB ecosystem"""
        ad_dtype = series.dtype
        print("================================================================ CMP", ad_dtype, tdb.dtype)
        # TileDB has no categorical, so assume it will convert to the type underlying the categorical
        if isinstance(ad_dtype, pd.CategoricalDtype):
            ad_dtype = series.cat.categories.dtype
        # TileDB has no object, so assume it will convert to the type underlying the object
        if ad_dtype == np.dtype("O"):
            ad_dtype = np.dtype(type(series[0]))
        # TileDB has no bool, and automatically converts to uint8
        if ad_dtype == bool:
            ad_dtype = np.uint8
        # xxxx foo
        if ad_dtype == np.float16:
            ad_dtype = np.float32

        print("---------------------------------------------------------------- NMP", ad_dtype, tdb.dtype)
        return ad_dtype == tdb.dtype

    for df_name in ["var", "obs"]:
        with tiledb.open((tmp_path / df_name).as_posix()) as arr:
            df = getattr(adata, df_name)
            # verify names match
            assert set(arr.schema.attr(i).name for i in range(arr.schema.nattr)) == set(
                getattr(adata, df_name).keys()
            )
            # verify length
            assert n == len(arr.query(dims=[]).df[:])
            # verify individual column types
            attr_idx = {
                arr.schema.attr(idx).name: idx for idx in range(arr.schema.nattr)
            }
            for k in df.keys():
                assert cmp_dtype(df[k], arr.schema.attr(attr_idx[k]))


# TODO: re-enable when #45 is resolved
@pytest.mark.skip(reason="Fails: filed as issue #45")
def test_from_anndata_annotations_empty(tmp_path):
    """
    Validate correct conversion with an empty (index-only) obs/var
    """
    n_obs = 100
    n_var = 10

    # AnnData requires a string index. TileDB does not support UTF8, so indices must be ASCII.
    obs = pd.DataFrame(index=np.arange(n_obs).astype(bytes))
    var = pd.DataFrame(index=np.arange(n_var).astype(bytes))

    X = np.ones((n_obs, n_var))
    adata = AnnData(X=X, obs=obs, var=var, dtype=X.dtype)

    SOMA(tmp_path.as_posix()).from_anndata(adata)

    assert all(
        (tmp_path / sub_array_path).exists()
        for sub_array_path in ["obs", "var", "X/data"]
    )

    # obs/var are sparse. Sort before comparing contents.
    with tiledb.open((tmp_path / "obs").as_posix()) as obs:
        assert np.array_equal(
            np.sort(adata.obs.index.to_numpy()), np.sort(obs[:]["obs_id"])
        )

    with tiledb.open((tmp_path / "var").as_posix()) as var:
        assert np.array_equal(
            np.sort(adata.var.index.to_numpy()), np.sort(var[:]["var_id"])
        )


# TODO: re-enable when #33 and #45 are resolved
@pytest.mark.skip(reason="Fails: filed as issues #33 and #45")
def test_from_anndata_annotations_none(tmp_path):
    """
    Validate ability to handle None in obs/var/X.
    """

    """ default constructor """
    path = tmp_path / "empty"
    adata = AnnData()
    SOMA(path.as_posix()).from_anndata(adata)
    assert all(
        (path / sub_array_path).exists() for sub_array_path in ["obs", "var", "X/data"]
    )

    """ only X defined """
    path = tmp_path / "X_only"
    adata = AnnData(X=np.eye(100, 10))
    SOMA(path.as_posix()).from_anndata(adata)
    assert all(
        (path / sub_array_path).exists() for sub_array_path in ["obs", "var", "X/data"]
    )

    """ missing var """
    path = tmp_path / "no_var"
    adata = AnnData(X=np.eye(100, 10), obs=np.arange(100))
    SOMA(path.as_posix()).from_anndata(adata)
    assert all(
        (path / sub_array_path).exists() for sub_array_path in ["obs", "var", "X/data"]
    )

    """ missing obs """
    path = tmp_path / "no_obs"
    adata = AnnData(X=np.eye(100, 10), var=np.arange(10))
    SOMA(path.as_posix()).from_anndata(adata)
    assert all(
        (path / sub_array_path).exists() for sub_array_path in ["obs", "var", "X/data"]
    )
