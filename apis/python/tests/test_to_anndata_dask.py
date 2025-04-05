from __future__ import annotations

from typing import Callable

import pytest

try:
    import dask.array as da

    DaskArray = da.Array
except ImportError:
    DaskArray = None

from scipy.sparse import csc_matrix, csr_matrix, hstack, vstack
from somacore import AxisQuery

from tiledbsoma import Experiment
from tiledbsoma._dask.load import load_daskarray
from tiledbsoma.io import to_anndata

from ._util import assert_array_equal, filter


@pytest.mark.parametrize("obs_chunk_size", [1, 5, 8, 20])
def test_dask_load_csr(
    conftest_pbmc_small_exp: Experiment,
    obs_chunk_size: int,
):
    layer = conftest_pbmc_small_exp.ms["RNA"].X["data"]
    X = load_daskarray(
        layer=layer,
        chunk_size=obs_chunk_size,
    )
    n_blocks = (X.shape[0] + obs_chunk_size - 1) // obs_chunk_size
    assert X.blocks.shape == (n_blocks, 1)
    blocks = [X.blocks[(i, 0)].compute() for i in range(n_blocks)]
    csr = vstack(blocks)
    X = X.compute()
    assert X.shape == (80, 20)
    assert X.nnz == 1600
    assert isinstance(X, csr_matrix)
    assert isinstance(csr, csr_matrix)
    assert_array_equal(X, csr)


@pytest.mark.parametrize("var_chunk_size", [1, 5, 8, 20])
def test_dask_load_csc(
    conftest_pbmc_small_exp: Experiment,
    var_chunk_size: int,
):
    layer = conftest_pbmc_small_exp.ms["RNA"].X["data"]
    X = load_daskarray(
        layer=layer,
        chunk_size=(None, var_chunk_size),
        format="csc",
    )
    n_blocks = (X.shape[1] + var_chunk_size - 1) // var_chunk_size
    assert X.blocks.shape == (1, n_blocks)
    blocks = [X.blocks[(0, i)].compute() for i in range(n_blocks)]
    csc = hstack(blocks)
    X = X.compute()
    assert X.shape == (80, 20)
    assert X.nnz == 1600
    assert isinstance(X, csc_matrix)
    assert isinstance(csc, csc_matrix)
    assert_array_equal(X, csc)


VerifyDaskArray = Callable[[csr_matrix, DaskArray], None]


@pytest.fixture
def verify_dask_array(
    obs_chunk_size: int,
    var_chunk_size: int,
    shape: tuple[int, int],
    nnz: int,
):
    """``fixture`` function that verifies a non-Dask AnnData matches a Dask-backed AnnData.

    Partially-applies other ``fixture``s, for calling convenience."""

    def check(
        X1: csr_matrix,
        X2: DaskArray,
    ):
        """Verify an AnnData (`ad1`) matches another (`ad2`) whose `X` matrix is a Dask Array."""
        assert X1.shape == shape
        assert X1.nnz == nnz
        nobs, nvar = X1.shape
        obs_chunks, var_chunks = X2.chunks
        assert len(obs_chunks) == (nobs + obs_chunk_size - 1) // obs_chunk_size
        assert len(var_chunks) == (nvar + var_chunk_size - 1) // var_chunk_size
        X2c = X2.compute()
        assert X2c.shape == shape
        assert X2c.nnz == nnz

        assert_array_equal(X1, X2c)

    return check


# fmt: off
sweep_queries = pytest.mark.parametrize(
    "obs_query,var_query,shape,nnz",
    [
        (AxisQuery(), AxisQuery(), (80, 20), 1600),
        (filter("nCount_RNA > 100"), AxisQuery(), (62, 20), 1240),
        (AxisQuery(), filter('attr("vst.variance.standardized") > 2.0'), (80, 5), 400),
        (filter("nCount_RNA > 100"), filter('attr("vst.variance.standardized") > 2.0'), (62, 5), 310),
    ],
)
# fmt: on


@sweep_queries
@pytest.mark.parametrize("obs_chunk_size", [20, 30, 80, 100])
@pytest.mark.parametrize("var_chunk_size", [3, 10, 20])
def test_dask_query_to_anndata(
    conftest_pbmc_small_exp: Experiment,
    obs_query: AxisQuery,
    var_query: AxisQuery,
    obs_chunk_size: int,
    var_chunk_size: int,
    verify_dask_array: VerifyDaskArray,
):
    query = conftest_pbmc_small_exp.axis_query(
        measurement_name="RNA",
        obs_query=obs_query,
        var_query=var_query,
    )
    ad1 = query.to_anndata(X_name="data")
    ad2 = query.to_anndata(
        X_name="data",
        dask=dict(chunk_size=(obs_chunk_size, var_chunk_size)),
    )
    verify_dask_array(ad1.X, ad2.X)


@pytest.mark.parametrize("obs_chunk_size", [20, 30, 80, 100])
@pytest.mark.parametrize("var_chunk_size", [3, 10, 20])
@pytest.mark.parametrize("shape", [(80, 20)])
@pytest.mark.parametrize("nnz", [1600])
def test_dask_experiment_to_anndata(
    conftest_pbmc_small_exp: Experiment,
    obs_chunk_size: int,
    var_chunk_size: int,
    verify_dask_array: VerifyDaskArray,
):
    ad1 = to_anndata(
        conftest_pbmc_small_exp,
        measurement_name="RNA",
    )
    ad2 = to_anndata(
        conftest_pbmc_small_exp,
        measurement_name="RNA",
        dask=dict(
            chunk_size=(obs_chunk_size, var_chunk_size),
        ),
    )
    verify_dask_array(ad1.X, ad2.X)


@sweep_queries
@pytest.mark.parametrize("obs_chunk_size", [20, 30, 80, 100])
@pytest.mark.parametrize("var_chunk_size", [3, 10, 20])
def test_sparse_read_dask_array(
    conftest_pbmc_small_exp: Experiment,
    obs_query: AxisQuery,
    var_query: AxisQuery,
    obs_chunk_size: int,
    var_chunk_size: int,
    verify_dask_array: VerifyDaskArray,
):
    query = conftest_pbmc_small_exp.axis_query(
        measurement_name="RNA",
        obs_query=obs_query,
        var_query=var_query,
    )
    ad1 = query.to_anndata(X_name="data")
    X2 = query.X("data").dask_array(chunk_size=(obs_chunk_size, var_chunk_size))
    verify_dask_array(ad1.X, X2)
