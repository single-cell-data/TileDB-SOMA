from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse
import tiledb

import tiledbsoma.io as io
from tiledbsoma import SOMA

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
        # "float16",    # TODO: Enable when #39 is fixed
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

    # AnnData requires string indices for obs/var
    obs = pd.DataFrame(
        data={"A": np.arange(n_obs, dtype=np.int32)}, index=np.arange(n_obs).astype(str)
    )
    var = pd.DataFrame(
        data={"A": np.arange(n_var, dtype=np.int32)}, index=np.arange(n_var).astype(str)
    )

    X_dtype = np.dtype(X_dtype_name)
    if X_encoding == "dense":
        X = np.eye(n_obs, n_var, dtype=X_dtype)
    elif X_encoding == "csc":
        X = scipy.sparse.eye(n_obs, n_var, dtype=X_dtype).tocsc()
    elif X_encoding == "csr":
        X = scipy.sparse.eye(n_obs, n_var, dtype=X_dtype).tocsr()

    adata = ad.AnnData(X=X, obs=obs, var=var, dtype=X.dtype)
    assert adata.X.dtype == X_dtype  # sanity

    io.from_anndata(SOMA(tmp_path.as_posix()), adata)
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
        (
            "categorical(int32)",
            lambda a: a.astype("int32").astype(
                pd.CategoricalDtype(categories=a.astype("int32"))
            ),
        ),
        (
            "categorical(uint32)",
            lambda a: a.astype("uint32").astype(
                pd.CategoricalDtype(categories=a.astype("uint32"))
            ),
        ),
        (
            "categorical(float32)",
            lambda a: a.astype("float32").astype(
                pd.CategoricalDtype(categories=a.astype("float32"))
            ),
        ),
        (
            "categorical(bool)",
            lambda a: a.astype("bool").astype(
                pd.CategoricalDtype(categories=a.astype("bool").unique())
            ),
        ),
    ]
    index = (
        np.arange(1, n + 1).astype(str).astype(bytes)
    )  # AnnData requires string indices, TileDB wants bytes. Use LCD
    df = pd.DataFrame(
        data={
            f"col_{name}": cast(pd.Series(index=index, data=np.arange(n)))
            for name, cast in df_col_type_sweep
        },
    )
    X = np.ones((n, n), dtype=np.float32)
    adata = ad.AnnData(X=X, obs=df, var=df, dtype=X.dtype)
    io.from_anndata(SOMA(tmp_path.as_posix()), adata)
    assert all(
        (tmp_path / sub_array_path).exists()
        for sub_array_path in ["obs", "var", "X/data"]
    )

    def cmp_dtype(series, tdb: tiledb.Attr) -> bool:
        """Encapsulate expected conversions moving into TileDB ecosystem"""
        ad_dtype = series.dtype
        # TileDB has no categorical, so assume it will convert to the type underlying the categorical
        if isinstance(ad_dtype, pd.CategoricalDtype):
            ad_dtype = series.cat.categories.dtype
        # TileDB has no object, so assume it will convert to the type underlying the object
        if ad_dtype == np.dtype("O"):
            ad_dtype = np.dtype(type(series[0]))
            # TODO: see annotation_dataframe.py. Once Unicode attributes are queryable, we'll need
            # to remove this check which is verifying the current force-to-ASCII workaround.
            if ad_dtype.name == "str":
                ad_dtype = np.dtype("S")

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

            # verify index
            assert np.array_equal(
                np.sort(df.index.to_numpy()), np.sort(arr[:][df_name + "_id"])
            )

            # verify individual column types
            attr_idx = {
                arr.schema.attr(idx).name: idx for idx in range(arr.schema.nattr)
            }
            for k in df.keys():
                assert cmp_dtype(df[k], arr.schema.attr(attr_idx[k]))


def test_from_anndata_annotations_empty(tmp_path):
    """
    Validate correct conversion with an empty (index-only) obs/var
    """
    n_obs = 100
    n_var = 10

    # AnnData requires a string index. TileDB does not support UTF8, so use ASCII.
    obs = pd.DataFrame(index=np.arange(n_obs).astype(bytes))
    var = pd.DataFrame(index=np.arange(n_var).astype(bytes))

    X = np.ones((n_obs, n_var))
    adata = ad.AnnData(X=X, obs=obs, var=var, dtype=X.dtype)

    io.from_anndata(SOMA(tmp_path.as_posix()), adata)

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


def test_from_anndata_annotations_none(tmp_path):
    """
    Validate correct handling of None in obs/var/X.
    """

    """ default constructor """
    path = tmp_path / "empty"
    adata = ad.AnnData()
    with pytest.raises(
        NotImplementedError, match="Empty AnnData.obs or AnnData.var unsupported."
    ):
        io.from_anndata(SOMA(path.as_posix()), adata)
    assert not any(
        (path / sub_array_path).exists() for sub_array_path in ["obs", "var", "X/data"]
    )

    """ only X defined """
    path = tmp_path / "X_only"
    adata = ad.AnnData(X=np.eye(100, 10, dtype=np.float32))
    io.from_anndata(SOMA(path.as_posix()), adata)
    assert all(
        (path / sub_array_path).exists() for sub_array_path in ["obs", "var", "X/data"]
    )

    """ missing var """
    path = tmp_path / "no_var"
    adata = ad.AnnData(
        X=np.eye(100, 10, dtype=np.float32), obs=np.arange(100).astype(str)
    )
    io.from_anndata(SOMA(path.as_posix()), adata)
    assert all(
        (path / sub_array_path).exists() for sub_array_path in ["obs", "var", "X/data"]
    )

    """ missing obs """
    path = tmp_path / "no_obs"
    adata = ad.AnnData(
        X=np.eye(100, 10, dtype=np.float32), var=np.arange(10).astype(str)
    )
    io.from_anndata(SOMA(path.as_posix()), adata)
    assert all(
        (path / sub_array_path).exists() for sub_array_path in ["obs", "var", "X/data"]
    )


def test_from_anndata_error_handling(tmp_path):
    """Ensure exception on a complex type we that should be unsupported."""
    n_obs = 10
    obs = pd.DataFrame(
        index=np.arange(n_obs).astype(str), data={"A": [{} for i in range(n_obs)]}
    )
    adata = ad.AnnData(obs=obs, X=np.ones((n_obs, 2), dtype=np.float32))
    with pytest.raises(NotImplementedError):
        io.from_anndata(SOMA(tmp_path.as_posix()), adata)


def test_from_anndata_zero_length_str(tmp_path):
    """
    Test case for issue #58: obs/var columns containing only zero length strings throw ArrowInvalid
    """
    n_obs = 100
    n_var = 10

    obs = pd.DataFrame(
        index=np.arange(n_obs).astype(str),
        data={
            "A": list(str(i) for i in range(n_obs)),
            "B": list("" for i in range(n_obs)),
        },
    )
    obs["A_cat"] = obs["A"].astype("category")
    obs["B_cat"] = obs["B"].astype("category")
    var = pd.DataFrame(
        index=np.arange(n_var).astype(str),
        data={"A": list(str(i) for i in range(n_var))},
    )
    X = np.ones((n_obs, n_var))
    adata = ad.AnnData(X=X, obs=obs, var=var, dtype=X.dtype)

    io.from_anndata(SOMA(tmp_path.as_posix()), adata)

    with tiledb.open((tmp_path / "obs").as_posix()) as obs:
        assert set(obs.schema.attr(i).name for i in range(obs.schema.nattr)) == set(
            adata.obs.keys()
        )
        assert adata.n_obs == len(obs.query(dims=[]).df[:])


test_nan_dtypes = [
    # (column_name, column_dtype, expect_raise)
    ("str", np.dtype(str), False),
    ("bytes", np.dtype(bytes), True),
    ("float64", np.float64, False),
    ("bool", np.dtype(bool), True),
    ("bool_", np.bool_, True),
    ("int64", np.int64, True),
    ("uint64", np.uint64, True),
]


@pytest.mark.parametrize("col_name,cat_dtype,expect_raise", test_nan_dtypes)
def test_from_anndata_category_nans(tmp_path, col_name, cat_dtype, expect_raise):
    """
    Categoricals can contain 'nan', ie, a series value which is a value not
    in the type's categories. While it conceptually represents "not a category",
    it is loosely referred to as a nan.

    Test conversion of various categorical series containing nans to ensure
    they are correctly handled.

    Presumed correct behavior depends on the underlying type of the category type,
    and follows the standard Pandas `astype()` coercion rules:
    * string: encode as 'nan'
    * float: encode as IEEE NaN
    * others: raise
    """
    n_obs = 8
    n_var = 4

    obs_idx = np.arange(n_obs).astype(str)
    obs = pd.DataFrame(
        index=obs_idx,
        data=pd.Categorical(
            np.arange(n_obs).astype(cat_dtype),
            categories=np.unique(np.arange(1, n_obs).astype(cat_dtype)),
        ),
        columns=[col_name],
    )
    var = pd.DataFrame(
        index=np.arange(n_var).astype(str),
        data=list(str(i) for i in range(n_var)),
        columns=["A"],
    )
    X = np.ones((n_obs, n_var), dtype=np.float32)
    adata = ad.AnnData(X=X, obs=obs, var=var)

    if expect_raise:
        with pytest.raises(ValueError):
            io.from_anndata(SOMA(tmp_path.as_posix()), adata)

    else:
        io.from_anndata(SOMA(tmp_path.as_posix()), adata)

        with tiledb.open((tmp_path / "obs").as_posix()) as arr:
            assert set(arr.schema.attr(i).name for i in range(arr.schema.nattr)) == set(
                adata.obs.keys()
            )
            obs_df = arr.df[:].sort_index()
            assert adata.n_obs == len(obs_df)
            assert np.array_equal(
                obs_df[col_name].astype(cat_dtype),
                adata.obs[col_name].astype(cat_dtype),
                equal_nan=True if np.dtype(cat_dtype).kind == "f" else False,
            )


def test_from_anndata_obsm_key_pandas_dataframe(tmp_path):
    """
    Normal case is:
    * X is scipy.sparse.csr_matrix or numpy.ndarray
    * obs,var are pandas.DataFrame (note _columns_ are numpy.ndarray)
    * obsm,varm elements are numpy.ndarray
    * obsp,varp elements are scipy.sparse.csr_matrix or numpy.ndarray
    Here we test the case where obsm has an element which is pandas.DataFrame

    See https://github.com/single-cell-data/TileDB-SingleCell/issues/74
    """
    input_path = Path(__file__).parent.parent / "anndata/pbmc-small.h5ad"
    adata = ad.read_h5ad(input_path)

    key = "is_a_pandas_dataframe"
    adata.obsm[key] = pd.DataFrame(data=np.zeros(adata.n_obs), index=adata.obs.index)
    adata.obsm[key].rename(columns={0: "column_name"}, inplace=True)

    io.from_anndata(SOMA(tmp_path.as_posix()), adata)
    assert all(
        (tmp_path / sub_array_path).exists()
        for sub_array_path in ["obs", "var", "X/data", "obsm", "obsm/" + key]
    )

    # Check types & shapes.
    adata_obsm = adata.obsm[key]
    with tiledb.open((tmp_path / "obsm" / key).as_posix()) as tiledb_obsm:
        assert adata_obsm.dtypes[0].kind == tiledb_obsm.schema.attr(0).dtype.kind
        assert adata_obsm.dtypes[0] == tiledb_obsm.schema.attr(0).dtype
        tiledb_obsm_df = tiledb_obsm[:]
        key_1 = list(tiledb_obsm_df.keys())[0]
        assert tiledb_obsm_df[key_1].shape == (adata.n_obs,)
