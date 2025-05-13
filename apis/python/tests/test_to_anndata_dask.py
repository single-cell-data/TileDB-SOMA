from __future__ import annotations

from pathlib import Path
from typing import Protocol

import numpy as np
import pytest

from .test_sparse_nd_array import create_random_tensor

try:
    import dask.array as da

    DaskArray = da.Array
except ImportError:
    DaskArray = None

import pyarrow as pa
from scipy.sparse import csc_matrix, csr_matrix, hstack, vstack
from somacore import AxisQuery

from tiledbsoma import Experiment, SparseNDArray
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


class VerifyDaskArray(Protocol):
    """Type-alias for the verifier-function returned by the ``verify_dask_array`` fixture below."""

    def __call__(
        self, X1: csr_matrix, X2: DaskArray, nnz: int | None = None
    ) -> None: ...


@pytest.fixture
def verify_dask_array(
    obs_chunk_size: int,
    var_chunk_size: int,
    shape: tuple[int, int],
    nnz: int,
) -> VerifyDaskArray:
    """``fixture`` function that verifies a non-Dask AnnData's ``X`` matches a Dask-backed AnnData's.

    Partially-applies other ``fixture``s, for calling convenience."""

    expected_nnz = nnz

    def check(
        X1: csr_matrix,
        X2: DaskArray,
        nnz: int | None = None,
    ):
        """Verify an AnnData (`ad1`) matches another (`ad2`) whose `X` matrix is a Dask Array."""
        if nnz is None:
            nnz = expected_nnz
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


@pytest.mark.parametrize(
    "obs_query,var_query,shape,nnz,obs_chunk_size,var_chunk_size",
    [(AxisQuery(), AxisQuery(), (80, 20), 1600, 20, 3)],
)
def test_dask_query_to_anndata_timestamp(
    conftest_pbmc_small_exp_path: Path,
    obs_query: AxisQuery,
    var_query: AxisQuery,
    obs_chunk_size: int,
    var_chunk_size: int,
    verify_dask_array: VerifyDaskArray,
):
    exp_path = str(conftest_pbmc_small_exp_path)

    def write_val(v: float) -> int:
        with Experiment.open(exp_path, mode="w") as exp:
            ts = exp.tiledb_timestamp_ms
            X = exp.ms["RNA"].X
            data = X["data"]
            shape = data.shape
            R, C = shape
            soma_dim_0, soma_dim_1, soma_data = zip(
                *[(r, c, v) for r in range(R) for c in range(C)]
            )
            tbl = pa.Table.from_pydict(
                dict(
                    soma_dim_0=soma_dim_0,
                    soma_dim_1=soma_dim_1,
                    soma_data=soma_data,
                )
            )
            data.write(tbl)
        return ts

    ts1 = write_val(111)
    ts2 = write_val(222)

    dask = dict(chunk_size=(obs_chunk_size, var_chunk_size))

    def check_ts(tiledb_timestamp: int | None, v: float):
        with Experiment.open(exp_path, tiledb_timestamp=tiledb_timestamp) as exp:
            data = exp.ms["RNA"].X["data"]
            [(data0, _)] = list(data.read().blockwise(0).scipy())
            data1 = data.read().dask_array(**dask).compute()
            assert_array_equal(data0, data1)
            assert data0[0, 0] == v

            query = exp.axis_query(
                measurement_name="RNA",
                obs_query=obs_query,
                var_query=var_query,
            )
            data = query.X("data")
            [(data0, _)] = data.blockwise(0).scipy()
            data1 = data.dask_array(**dask).compute()
            assert_array_equal(data0, data1)
            assert data0[0, 0] == v

            ad1 = to_anndata(exp, measurement_name="RNA")
            ad2 = to_anndata(exp, measurement_name="RNA", dask=dask)
            verify_dask_array(ad1.X, ad2.X)

            ad1 = query.to_anndata(X_name="data")
            ad2 = query.to_anndata(X_name="data", dask=dask)
            verify_dask_array(ad1.X, ad2.X)

    check_ts(ts1 - 1, -0.1934951510503384)
    check_ts(ts1, 111)
    check_ts(ts2 - 1, 111)
    check_ts(ts2, 222)
    check_ts(None, 222)


@pytest.mark.parametrize(
    "obs_query,var_query,shape,nnz,obs_chunk_size,var_chunk_size",
    [(AxisQuery(), AxisQuery(), (80, 20), 1600, 20, 3)],
)
def test_dask_query_to_anndata_layers(
    conftest_pbmc_small_exp_path: Path,
    obs_query: AxisQuery,
    var_query: AxisQuery,
    obs_chunk_size: int,
    var_chunk_size: int,
    verify_dask_array: VerifyDaskArray,
):
    exp_path = str(conftest_pbmc_small_exp_path)
    layer_nnz = 528
    with Experiment.open(exp_path, mode="w") as exp:
        X = exp.ms["RNA"].X
        shape = X["data"].shape
        layer1 = SparseNDArray.create(
            f"{X.uri}/layer1",
            type=pa.float32(),
            shape=shape,
            context=X.context,
        )
        arrow_tensor = create_random_tensor(
            format="csr",
            shape=shape,
            dtype=np.float32(),
            seed=1111,
        )
        assert arrow_tensor.non_zero_length == layer_nnz
        layer1.write(arrow_tensor)
        X.set("layer1", layer1)

    with Experiment.open(exp_path) as exp:
        query = exp.axis_query(
            measurement_name="RNA",
            obs_query=obs_query,
            var_query=var_query,
        )
        ad1 = query.to_anndata(X_name="data", X_layers=("layer1",))
        ad2 = query.to_anndata(
            X_name="data",
            X_layers=("layer1",),
            dask=dict(chunk_size=(obs_chunk_size, var_chunk_size)),
        )
        verify_dask_array(ad1.X, ad2.X)
        verify_dask_array(ad1.layers["layer1"], ad2.layers["layer1"], nnz=layer_nnz)


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
