import numpy as np
import pyarrow as pa
import pytest
from scipy import sparse as sp

import tiledbsoma as soma
import tiledbsoma.io as somaio


@pytest.fixture
def src_matrix(request):
    format, shape, density = request.param
    if format == "dense":
        return (
            np.random.default_rng()
            .standard_normal(np.prod(shape), dtype=np.float32)
            .reshape(shape)
        )

    return sp.random(10, 89, density=density, format=format, dtype=np.float32)


@pytest.mark.parametrize(
    "plfm_config",
    [
        soma.TileDBPlatformConfig(write_X_chunked=False, goal_chunk_nnz=10000),
        soma.TileDBPlatformConfig(write_X_chunked=False, goal_chunk_nnz=100000),
        soma.TileDBPlatformConfig(write_X_chunked=True, goal_chunk_nnz=10000),
        soma.TileDBPlatformConfig(write_X_chunked=True, goal_chunk_nnz=100000),
    ],
)
@pytest.mark.parametrize(
    "src_matrix",
    [
        ("dense", (10, 100), 1),
        ("dense", (1103, 107), 1),
        ("csc", (10, 89), 0.01),
        ("csr", (10, 89), 0.01),
        ("csc", (1001, 899), 0.3),
        ("csr", (1001, 899), 0.3),
    ],
    indirect=True,
)
def test_io_create_from_matrix_Dense_nd_array(tmp_path, plfm_config, src_matrix):
    """
    Test soma.io.from_matrix to a DenseNdArray.

    Cases:
    * src matrix is:  csc, csr, ndarray
    * _tiledb_platform_config.write_X_chunked: True or False
    * src_array bigger or smaller than _tiledb_platform_config.goal_chunk_nnz
    """
    snda = soma.DenseNdArray(
        tmp_path.as_posix(), tiledb_platform_config=plfm_config
    )
    somaio.create_from_matrix(snda, src_matrix)

    assert snda.shape == src_matrix.shape
    assert snda.ndim == src_matrix.ndim

    read_back = snda.read_numpy((slice(None), slice(None)))

    if isinstance(src_matrix, np.ndarray):
        assert np.array_equal(read_back, src_matrix)
    else:
        assert np.array_equal(read_back, src_matrix.toarray())


@pytest.mark.parametrize(
    "plfm_config",
    [
        soma.TileDBPlatformConfig(write_X_chunked=False, goal_chunk_nnz=10000),
        soma.TileDBPlatformConfig(write_X_chunked=False, goal_chunk_nnz=100000),
        soma.TileDBPlatformConfig(write_X_chunked=True, goal_chunk_nnz=10000),
        soma.TileDBPlatformConfig(write_X_chunked=True, goal_chunk_nnz=100000),
    ],
)
@pytest.mark.parametrize(
    "src_matrix",
    [
        ("dense", (10, 100), 1),
        ("dense", (1103, 107), 1),
        ("csc", (10, 89), 0.01),
        ("csr", (10, 89), 0.01),
        ("csc", (1001, 899), 0.3),
        ("csr", (1001, 899), 0.3),
    ],
    indirect=True,
)
def test_io_create_from_matrix_Sparse_nd_array(tmp_path, plfm_config, src_matrix):
    """
    Test soma.io.from_matrix to a SparseNdArray.

    Cases:
    * src matrix is:  csc, csr, ndarray
    * _tiledb_platform_config.write_X_chunked: True or False
    * src_array bigger or smaller than _tiledb_platform_config.goal_chunk_nnz
    """
    snda = soma.SparseNdArray(
        tmp_path.as_posix(), tiledb_platform_config=plfm_config
    )
    somaio.create_from_matrix(snda, src_matrix)

    assert snda.shape == src_matrix.shape
    assert snda.ndim == src_matrix.ndim

    tbl = pa.concat_tables(snda.read_table((slice(None), slice(None))))
    read_back = sp.csr_matrix(
        (
            tbl.column("soma_data").to_numpy(),
            (tbl.column("soma_dim_0").to_numpy(), tbl.column("soma_dim_1").to_numpy()),
        ),
        shape=snda.shape,
    )

    # fast equality check using __ne__
    assert (sp.csr_matrix(src_matrix) != read_back).nnz == 0
