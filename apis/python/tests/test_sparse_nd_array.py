from __future__ import annotations

import operator
import pathlib
import sys
from typing import Any, Dict, List, Tuple, Union

import numpy as np
import pyarrow as pa
import pytest
import scipy.sparse as sparse
import tiledb

import tiledbsoma as soma
from tiledbsoma import _factory
from tiledbsoma.options import SOMATileDBContext

from . import NDARRAY_ARROW_TYPES_NOT_SUPPORTED, NDARRAY_ARROW_TYPES_SUPPORTED

AnySparseTensor = Union[pa.SparseCOOTensor, pa.SparseCSRMatrix, pa.SparseCSCMatrix]


@pytest.mark.parametrize(
    "shape", [(10,), (1, 100), (10, 1, 100), (2, 4, 6, 8), [1], (1, 2, 3, 4, 5)]
)
@pytest.mark.parametrize("element_type", NDARRAY_ARROW_TYPES_SUPPORTED)
def test_sparse_nd_array_create_ok(
    tmp_path, shape: Tuple[int, ...], element_type: pa.DataType
):
    """
    Test all cases we expect "create" to succeed.
    """
    assert pa.types.is_primitive(element_type)  # sanity check incoming params

    with pytest.raises(TypeError):
        # non-arrow write
        soma.SparseNDArray.create(
            tmp_path.as_posix(), type=element_type.to_pandas_dtype(), shape=shape
        )
    a = soma.SparseNDArray.create(tmp_path.as_posix(), type=element_type, shape=shape)
    assert soma.SparseNDArray.exists(tmp_path.as_posix())
    assert not soma.DenseNDArray.exists(tmp_path.as_posix())
    assert not soma.Collection.exists(tmp_path.as_posix())
    assert a.soma_type == "SOMASparseNDArray"
    assert a.uri == tmp_path.as_posix()
    assert a.ndim == len(shape)
    assert a.shape == tuple(shape)
    assert a.is_sparse is True

    assert a.schema is not None
    expected_field_names = ["soma_data"] + [f"soma_dim_{d}" for d in range(len(shape))]
    assert set(a.schema.names) == set(expected_field_names)
    for d in range(len(shape)):
        assert a.schema.field(f"soma_dim_{d}").type == pa.int64()
    assert a.schema.field("soma_data").type == element_type


@pytest.mark.parametrize("shape", [(10,)])
@pytest.mark.parametrize("element_type", NDARRAY_ARROW_TYPES_NOT_SUPPORTED)
def test_sparse_nd_array_create_fail(
    tmp_path, shape: Tuple[int, ...], element_type: pa.DataType
):
    with pytest.raises(TypeError):
        soma.SparseNDArray.create(tmp_path.as_posix(), type=element_type, shape=shape)


def create_random_tensor(
    format: str,
    shape: Tuple[int, ...],
    dtype: np.dtype,
    density: float = 0.33,
):
    """
    Create a random tensor/table of specified format, shape and dtype.

    Guarantee: there will be NO duplicate values in the tensor, which makes
    it simpler to validate (see `tensors_are_same_value`)
    """
    rng = np.random.default_rng()
    ndim = len(shape)
    assert 0 < density <= 1

    assert format in ["coo", "csc", "csr", "table"], "Unimplemented format"

    if format == "coo":
        nrec = int(density * np.prod(shape))
        data = rng.choice(10 * nrec, size=nrec, replace=False).astype(dtype)
        all_coords = np.array(
            np.meshgrid(*tuple(np.arange(dim_len) for dim_len in shape))
        ).T.reshape(-1, ndim)
        coords = rng.choice(all_coords, nrec, replace=False)
        return pa.SparseCOOTensor.from_numpy(data, coords, shape=shape)

    if format == "table":
        nrec = int(density * np.prod(shape))
        data = rng.choice(10 * nrec, size=nrec, replace=False).astype(dtype)
        all_coords = np.array(
            np.meshgrid(*tuple(np.arange(dim_len) for dim_len in shape))
        ).T.reshape(-1, ndim)
        coords = rng.choice(all_coords, nrec, replace=False).T
        pydict = {
            f"soma_dim_{n}": pa.array(coords[n], type=pa.int64()) for n in range(ndim)
        }
        pydict.update({"soma_data": pa.array(data)})
        return pa.Table.from_pydict(pydict)

    if format == "csc":
        assert ndim == 2
        return pa.SparseCSCMatrix.from_scipy(
            sparse.random(
                *shape,
                density=density,
                format=format,
                random_state=rng,
                dtype=dtype,
            )
        )

    if format == "csr":
        assert ndim == 2
        return pa.SparseCSRMatrix.from_scipy(
            sparse.random(
                *shape,
                density=density,
                format=format,
                random_state=rng,
                dtype=dtype,
            )
        )


def tensors_are_same_value(a: AnySparseTensor, b: AnySparseTensor) -> bool:
    """
    Return true if the tenors contain the same values, allowing for
    differences in coordinate ordering
    """
    if type(a) != type(b):
        return False
    if a.shape != b.shape:
        return False
    if a.type != b.type:
        return False

    def _check_coo_values(a, b) -> bool:
        a_data, a_coords = a
        b_data, b_coords = b
        ai = a_data.flatten().argsort().reshape(-1, 1)
        bi = b_data.flatten().argsort().reshape(-1, 1)

        if not np.array_equal(
            np.take_along_axis(a_data, ai, axis=0),
            np.take_along_axis(b_data, bi, axis=0),
        ):
            return False
        if not np.array_equal(
            np.take_along_axis(a_coords, ai, axis=0),
            np.take_along_axis(b_coords, bi, axis=0),
        ):
            return False
        return True

    # coordinate order in the tensors may not be the same, leading to these gymnastics
    if isinstance(a, pa.SparseCOOTensor):
        return _check_coo_values(a.to_numpy(), b.to_numpy())

    if isinstance(a, (pa.SparseCSRMatrix, pa.SparseCSCMatrix)):
        _a = pa.SparseCOOTensor.from_scipy(a.to_scipy().tocoo())
        _b = pa.SparseCOOTensor.from_scipy(b.to_scipy().tocoo())
        return _check_coo_values(_a.to_numpy(), _b.to_numpy())

    return False


def tables_are_same_value(a: pa.Table, b: pa.Table) -> bool:
    """
    Return True if the tables contain the same COO array data,
    allowing for differences in coordinate order.
    """
    if a.shape != b.shape:
        return False
    if a.field("soma_data").type != b.field("soma_data").type:
        return False
    for tbl in (a, b):
        if not all(
            tbl.field(f"soma_dim_{n}").type == pa.int64()
            for n in range(tbl.num_columns - 1)
        ):
            return False

    ndim = a.shape[1] - 1
    ai = a.column("soma_data").to_numpy().argsort()
    bi = b.column("soma_data").to_numpy().argsort()
    if not np.array_equal(
        np.take_along_axis(a.column("soma_data").to_numpy(), ai, axis=0),
        np.take_along_axis(b.column("soma_data").to_numpy(), bi, axis=0),
    ):
        return False

    for n in range(ndim):
        dim_name = f"soma_dim_{n}"
        if a.field(dim_name).type != pa.int64() or b.field(dim_name).type != pa.int64():
            return False
        if not np.array_equal(
            np.take_along_axis(a.column(dim_name).to_numpy(), ai, axis=0),
            np.take_along_axis(b.column(dim_name).to_numpy(), bi, axis=0),
        ):
            return False

    return True


@pytest.mark.parametrize(
    "shape,format",
    [
        ((10,), "coo"),
        ((10, 21), "coo"),
        ((10, 20, 2), "coo"),
        ((2, 4, 6, 8), "coo"),
        ((1, 2, 4, 6, 8), "coo"),
    ],
)
@pytest.mark.parametrize("test_enumeration", range(10))
def test_sparse_nd_array_read_write_sparse_tensor(
    tmp_path,
    shape: Tuple[int, ...],
    format: str,
    test_enumeration: int,
):
    # Test sanity: Tensor only, and CSC and CSR only support 2D, so fail any nonsense configs
    assert format == "coo"

    a = soma.SparseNDArray.create(tmp_path.as_posix(), type=pa.float64(), shape=shape)
    assert a.shape == shape

    # Make a random sample in the desired format
    # As discussed in the SparseNDArray implementation, Arrow SparseTensor objects can't be zero-length
    # so we must be prepared for StopIteration on reading them. It simplifies unit-test logic to use
    # occupation density of 1.0 for this test.
    data = create_random_tensor(format, shape, np.float64, 1.0)
    with pytest.raises(TypeError):
        # non-arrow write
        a.write(data.to_numpy())
    a.write(data)
    del a

    # Read back and validate
    with soma.SparseNDArray.open(tmp_path.as_posix()) as b:
        t = b.read((slice(None),) * len(shape)).coos().concat()

        assert tensors_are_same_value(t, data)

        t = next(b.read((0,) * len(shape)).coos())

        assert t.shape == shape

    # Validate TileDB array schema
    with tiledb.open(tmp_path.as_posix()) as A:
        assert A.schema.sparse
        assert not A.schema.allows_duplicates


@pytest.mark.parametrize("shape", [(10,), (23, 4), (5, 3, 1), (8, 4, 2, 30)])
@pytest.mark.parametrize("test_enumeration", range(10))
def test_sparse_nd_array_read_write_table(
    tmp_path, shape: Tuple[int, ...], test_enumeration: int
):
    a = soma.SparseNDArray.create(tmp_path.as_posix(), type=pa.float32(), shape=shape)
    assert a.shape == shape

    # make a random sample in the desired format
    data = create_random_tensor("table", shape, np.float32)
    a.write(data)
    del a

    # Read back and validate
    b = soma.SparseNDArray.open(tmp_path.as_posix())
    t = next(b.read((slice(None),) * len(shape)).tables())
    assert isinstance(t, pa.Table)
    assert tables_are_same_value(data, t)

    # Validate TileDB array schema
    with tiledb.open(tmp_path.as_posix()) as A:
        assert A.schema.sparse
        assert not A.schema.allows_duplicates


@pytest.mark.parametrize("dtype", [np.float32, np.float64, np.int32, np.int64])
@pytest.mark.parametrize("shape", [(1,), (23, 14), (35, 3, 2), (8, 4, 2, 30)])
def test_sparse_nd_array_read_as_pandas(
    tmp_path, dtype: np.dtype, shape: Tuple[int, ...]
):
    dtype = np.dtype(dtype)
    with soma.SparseNDArray.create(
        tmp_path.as_posix(), type=pa.from_numpy_dtype(dtype), shape=shape
    ) as a:
        assert a.shape == shape

        # make a random sample in the desired format
        data = create_random_tensor("table", shape, dtype)
        a.write(data)

    with _factory.open(tmp_path.as_posix()) as a:
        df = a.read().tables().concat().to_pandas()

    dim_names = [f"soma_dim_{n}" for n in range(len(shape))]
    assert df.sort_values(by=dim_names, ignore_index=True).equals(
        data.to_pandas().sort_values(by=dim_names, ignore_index=True)
    )

    # Validate TileDB array schema
    with tiledb.open(tmp_path.as_posix()) as A:
        assert A.schema.sparse
        assert not A.schema.allows_duplicates


@pytest.mark.parametrize("shape_is_nones", [True, False])
@pytest.mark.parametrize("element_type", NDARRAY_ARROW_TYPES_SUPPORTED)
def test_sparse_nd_array_shaping(tmp_path, shape_is_nones, element_type):
    uri = tmp_path.as_posix()

    soma.SparseNDArray.create(
        uri,
        type=element_type,
        shape=(None, None) if shape_is_nones else (2, 3),
    ).close()
    assert soma.SparseNDArray.exists(uri)

    with soma.SparseNDArray.open(uri) as snda:
        assert snda.nnz == 0

    batch1 = pa.Table.from_pydict(
        {
            "soma_dim_0": [0, 0, 0, 1, 1, 1],
            "soma_dim_1": [0, 1, 2, 0, 1, 2],
            "soma_data": [1, 2, 3, 4, 5, 6],
        }
    )

    batch2 = pa.Table.from_pydict(
        {"soma_dim_0": [2, 2, 2], "soma_dim_1": [0, 1, 2], "soma_data": [7, 8, 9]}
    )

    with soma.SparseNDArray.open(uri, "w") as snda:
        snda.write(batch1)

    with soma.SparseNDArray.open(uri) as snda:
        assert snda.nnz == 6

    if shape_is_nones:
        with soma.SparseNDArray.open(uri, "w") as snda:
            snda.write(batch2)
    else:
        # tiledb.cc.TileDBError: [TileDB::Dimension] Error: Coordinate 2
        # is out of domain bounds [0, 1] on dimension 'soma_dim_0'
        with pytest.raises(tiledb.cc.TileDBError):
            with soma.SparseNDArray.open(uri, "w") as snda:
                snda.write(batch2)


def test_empty_read(tmp_path):
    """
    Verify that queries expected to return empty results actually
    work. There are edge cases around SparseTensors, which are unable
    to represent empty arrays.
    """
    soma.SparseNDArray.create(
        tmp_path.as_posix(), type=pa.uint16(), shape=(10, 100)
    ).close()

    with soma.SparseNDArray.open(tmp_path.as_posix()) as a:
        #
        # First, test reads of zero element array
        #

        # These work as expected
        coords = (slice(None),)
        assert sum(len(t) for t in a.read(coords).tables()) == 0
        # Fails due to ARROW-17933
        # assert sum(t.non_zero_length for t in a.read(coords).coos()) == 0

    #
    # Next, test empty queries on non-empty array
    #
    with soma.SparseNDArray.open(tmp_path.as_posix(), "w") as a:
        a.write(
            pa.SparseCOOTensor.from_scipy(
                sparse.coo_matrix(([1], ([0], [0])), shape=a.shape)
            )
        )
    with soma.SparseNDArray.open(tmp_path.as_posix()) as a:
        assert sum(len(t) for t in a.read((slice(None),)).tables()) == 1

        coords = (1, 1)  # no element at this coordinate
        assert sum(len(t) for t in a.read(coords).tables()) == 0


@pytest.mark.xfail(sys.version_info.minor > 7, reason="bug ARROW-17933")
def test_empty_read_sparse_coo(tmp_path):
    """
    this test is factored from test_empty_read() because it is subject
    to ARROW-17933, and is xfailed on certain python verisons. The tests
    can be recombined with test_empty_read once the behavior is consistent
    across versions.

    ---

    TODO: Due to bug https://issues.apache.org/jira/browse/ARROW-17933, this
    _incorrectly_ raises an ArrowInvalid exception. Remove the `pyarrow.throws`
    when fixed, as it is supported API and should work.

    It does NOT fail on Python3.7, but does fail on later versions (unclear why,
    perhaps a NumPy difference)

    """
    soma.SparseNDArray.create(
        tmp_path.as_posix(), type=pa.uint16(), shape=(10, 100)
    ).close()

    with soma.SparseNDArray.open(tmp_path.as_posix()) as a:
        coords = (slice(None),)
        assert sum(t.non_zero_length for t in a.read(coords).coos()) == 0

    with soma.SparseNDArray.open(tmp_path.as_posix(), "w") as a:
        a.write(
            pa.SparseCOOTensor.from_scipy(
                sparse.coo_matrix(([1], ([0], [0])), shape=a.shape)
            )
        )
    with soma.SparseNDArray.open(tmp_path.as_posix()) as a:
        assert sum(len(t) for t in a.read((slice(None),)).tables()) == 1

        coords = (1, 1)  # no element at this coordinate
        assert sum(t.non_zero_length for t in a.read(coords).coos()) == 0


def test_coo_custom_shape(tmp_path):
    soma.SparseNDArray.create(
        tmp_path.as_posix(), type=pa.uint16(), shape=(1000, 1000)
    ).close()

    with soma.SparseNDArray.open(tmp_path.as_posix()) as a:
        coords = (slice(None),)
        assert a.read(coords).coos().shape == (1000, 1000)
        assert a.read(coords).coos(shape=(500, 500)).shape == (500, 500)
        with pytest.raises(ValueError):
            assert a.read(coords).coos(shape=(500, 500, 500))


@pytest.mark.parametrize("shape", [(), (0,), (10, 0), (0, 10), (1, 2, 0)])
def test_zero_length_fail(tmp_path, shape):
    """Zero length dimensions are expected to fail"""
    with pytest.raises(ValueError):
        soma.SparseNDArray.create(tmp_path.as_posix(), type=pa.float32(), shape=shape)


def test_sparse_nd_array_nnz(tmp_path):
    soma.SparseNDArray.create(
        tmp_path.as_posix(), type=pa.int32(), shape=(10, 10, 10)
    ).close()
    with soma.SparseNDArray.open(tmp_path.as_posix()) as a:
        assert a.nnz == 0
    with soma.SparseNDArray.open(tmp_path.as_posix(), "w") as a:
        t: pa.SparseCOOTensor = create_random_tensor(
            "coo", a.shape, pa.int32().to_pandas_dtype(), 0.1
        )
        a.write(t)
    with soma.SparseNDArray.open(tmp_path.as_posix()) as a:
        assert t.non_zero_length == a.nnz


def test_sparse_nd_array_reshape(tmp_path):
    """
    Reshape currently unimplemented.
    """
    with soma.SparseNDArray.create(
        tmp_path.as_posix(), type=pa.int32(), shape=(10, 10, 10)
    ) as a:
        with pytest.raises(NotImplementedError):
            assert a.reshape((100, 10, 1))


@pytest.mark.parametrize(
    "shape",
    [(4,), (4, 5, 6)],
)
def test_csr_csc_2d_read(tmp_path, shape):
    """Arrays which are not 2D can't be requested in CSC or CSR format."""

    arrow_tensor = create_random_tensor(
        format="coo",
        shape=shape,
        dtype=np.float32(),
    )

    snda = soma.SparseNDArray.create(
        tmp_path.as_posix(), type=pa.float64(), shape=shape
    )
    snda.write(arrow_tensor)


@pytest.mark.parametrize(
    "write_format",
    ["coo", "csr", "csc"],
)
@pytest.mark.parametrize(
    # We want to test read_format == "none_of_the_above", to ensure it throws NotImplementedError,
    # but that can't be gotten past typeguard.
    "read_format",
    ["table", "coo"],
)
@pytest.mark.parametrize(
    "io",
    [
        {
            "name": "coords=()",
            "shape": (4,),
            "coords": (),
            "dims": {
                "soma_dim_0": [0, 1, 2, 3],
            },
            "throws": None,
        },
        {
            "name": "coords=[None]",
            "shape": (4,),
            "coords": (None,),
            "dims": {
                "soma_dim_0": [0, 1, 2, 3],
            },
            "throws": None,
        },
        {
            "name": "coords=[[-100:100]]",
            "shape": (4,),
            "coords": (slice(-100, 100),),
            "dims": {
                "soma_dim_0": [0, 1, 2, 3],
            },
            "throws": None,
        },
        {
            "name": "coords=[4]",
            "shape": (4,),
            "coords": (1,),
            "dims": {
                "soma_dim_0": [1],
            },
            "throws": None,
        },
        {
            "name": "coords=[[2, 4]]",
            "shape": (6,),
            "coords": [[2, 4]],
            "dims": {
                "soma_dim_0": [2, 4],
            },
            "throws": None,
        },
        {
            "name": "negative unsupported",
            "shape": (6,),
            "coords": [[-2, -4]],
            "dims": {
                "soma_dim_0": [2, 4],
            },
            "throws": (
                RuntimeError,
                tiledb.cc.TileDBError,
            ),
        },
        {
            "name": "coords=[0,0]",
            "shape": (4, 6),
            "coords": (0, 0),
            "dims": {
                "soma_dim_0": [0],
                "soma_dim_1": [0],
            },
            "throws": None,
        },
        {
            "name": "coords=([:2],[2:])",
            "shape": (3, 4),
            "coords": (slice(None, 2), slice(2, None)),
            "dims": {
                "soma_dim_0": [0, 0, 1, 1, 2, 2],
                "soma_dim_1": [2, 3, 2, 3, 2, 3],
            },
            "throws": None,
        },
        {
            "name": "coords=([-10:2],[2:500])",
            "shape": (3, 4),
            "coords": (slice(-10, 2), slice(2, 500)),
            "dims": {
                "soma_dim_0": [0, 0, 1, 1, 2, 2],
                "soma_dim_1": [2, 3, 2, 3, 2, 3],
            },
            "throws": None,
        },
        {
            "name": "2D coords=[0]",
            "shape": (4, 6),
            "coords": (0,),  # Remaining dimensions are implicit-all
            "dims": {
                "soma_dim_0": [0, 0, 0, 0, 0, 0],
                "soma_dim_1": [0, 1, 2, 3, 4, 5],
            },
            "throws": None,
        },
        {
            "name": "too many coords",
            "shape": (4, 6),
            "coords": (0, 0, 0),
            "dims": {
                "soma_dim_0": [0, 0, 0, 0, 0, 0],
                "soma_dim_1": [0, 1, 2, 3, 4, 5],
            },
            "throws": ValueError,
        },
        {
            "name": "3D coords",
            "shape": (4, 5, 6),
            "coords": (2, 3, 4),
            "dims": {
                "soma_dim_0": [2],
                "soma_dim_1": [3],
                "soma_dim_2": [4],
            },
            "throws": None,
        },
        {
            "name": "2D coords=(3, 4)",
            "shape": (4, 6),
            "coords": (3, 4),
            "dims": {
                "soma_dim_0": [3],
                "soma_dim_1": [4],
            },
            "throws": None,
        },
        {
            "name": "2D coords=([1:2], [3:4])",
            "shape": (4, 6),
            "coords": (slice(1, 2), slice(3, 4)),
            "dims": {
                "soma_dim_0": [1, 1, 2, 2],
                "soma_dim_1": [3, 4, 3, 4],
            },
            "throws": None,
        },
        {
            "name": "2D coords=([1:2], [3, 4])",
            "shape": (4, 6),
            "coords": (slice(1, 2), [3, 4]),
            "dims": {
                "soma_dim_0": [1, 1, 2, 2],
                "soma_dim_1": [3, 4, 3, 4],
            },
            "throws": None,
        },
        {
            "name": "2D coords=(np[1, 2], pa[3, 4])",
            "shape": (4, 6),
            "coords": (np.asarray([1, 2]), pa.array([3, 4])),
            "dims": {
                "soma_dim_0": [1, 1, 2, 2],
                "soma_dim_1": [3, 4, 3, 4],
            },
            "throws": None,
        },
        {
            "name": "2D coords=(np[[1, 2]], pa[3, 4])",
            "shape": (4, 6),
            "coords": (np.asarray([[1, 2]]), pa.array([3, 4])),
            "dims": {
                "soma_dim_0": [1, 1, 2, 2],
                "soma_dim_1": [3, 4, 3, 4],
            },
            "throws": ValueError,  # np.ndarray must be 1D
        },
        {
            "name": "2D coords=([:], [3:4])",
            "shape": (4, 6),
            "coords": (slice(None), slice(3, 4)),
            "dims": {
                "soma_dim_0": [0, 0, 1, 1, 2, 2, 3, 3],
                "soma_dim_1": [3, 4, 3, 4, 3, 4, 3, 4],
            },
            "throws": None,
        },
        {
            "name": "2D coords=([:], [3:100])",
            "shape": (4, 6),
            "coords": (slice(None), slice(3, 100)),
            "dims": {
                "soma_dim_0": [0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3],
                "soma_dim_1": [3, 4, 5, 3, 4, 5, 3, 4, 5, 3, 4, 5],
            },
            "throws": None,
        },
        {
            "name": "2D coords=([1:2], [:])",
            "shape": (4, 6),
            "coords": (slice(1, 2), slice(None)),
            "dims": {
                "soma_dim_0": [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2],
                "soma_dim_1": [0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5],
            },
            "throws": None,
        },
        {
            "name": "2D coords=([:], [:])",
            "shape": (3, 4),
            "coords": (slice(None), slice(None)),
            "dims": {
                "soma_dim_0": [
                    0,
                    0,
                    0,
                    0,
                    1,
                    1,
                    1,
                    1,
                    2,
                    2,
                    2,
                    2,
                ],
                "soma_dim_1": [
                    0,
                    1,
                    2,
                    3,
                    0,
                    1,
                    2,
                    3,
                    0,
                    1,
                    2,
                    3,
                ],
            },
            "throws": None,
        },
        {
            "name": "2D coords=[]",
            "shape": (3, 4),
            "coords": [],
            "dims": {
                "soma_dim_0": [
                    0,
                    0,
                    0,
                    0,
                    1,
                    1,
                    1,
                    1,
                    2,
                    2,
                    2,
                    2,
                ],
                "soma_dim_1": [
                    0,
                    1,
                    2,
                    3,
                    0,
                    1,
                    2,
                    3,
                    0,
                    1,
                    2,
                    3,
                ],
            },
            "throws": None,
        },
        {
            "name": "3D coords=([1:2], [2:3], [3:4])",
            "shape": (4, 5, 6),
            "coords": (slice(1, 2), slice(2, 3), slice(3, 4)),
            "dims": {
                "soma_dim_0": [1, 1, 1, 1, 2, 2, 2, 2],
                "soma_dim_1": [2, 2, 3, 3, 2, 2, 3, 3],
                "soma_dim_2": [3, 4, 3, 4, 3, 4, 3, 4],
            },
            "throws": None,
        },
        {
            "name": "2D coords=(np32[1, 2], np64[3, 4])",
            "shape": (9, 11),
            "coords": (
                np.array([1, 2], dtype=np.int32),
                np.array([3, 4], dtype=np.int64),
            ),
            "dims": {
                "soma_dim_0": [1, 1, 2, 2],
                "soma_dim_1": [3, 4, 3, 4],
            },
            "throws": None,
        },
        {
            "name": "2D coords=(np64[1, 2], np32[3, 4])",
            "shape": (9, 11),
            "coords": (
                np.array([1, 2], dtype=np.uint64),
                np.array([3, 4], dtype=np.uint32),
            ),
            "dims": {
                "soma_dim_0": [1, 1, 2, 2],
                "soma_dim_1": [3, 4, 3, 4],
            },
            "throws": None,
        },
    ],
    ids=lambda io: io.get("name"),
)
def test_sparse_nd_array_table_slicing(tmp_path, io, write_format, read_format):
    if (write_format == "csr" or write_format == "csc") and len(io["shape"]) != 2:
        return  # Not supported by create_random_tensor
    if (read_format == "csr" or read_format == "csc") and len(io["shape"]) != 2:
        return  # Not supported by readback; exception-throwing for this is tested separately above.

    # Set up contents
    arrow_tensor = create_random_tensor(
        format=write_format,
        shape=io["shape"],
        dtype=np.float32(),
        density=1.0,
    )

    with soma.SparseNDArray.create(
        tmp_path.as_posix(), type=pa.float64(), shape=io["shape"]
    ) as snda_w:
        snda_w.write(arrow_tensor)

    with soma.SparseNDArray.open(tmp_path.as_posix()) as snda:
        if read_format == "table":
            if io["throws"] is not None:
                with pytest.raises(io["throws"]):
                    next(snda.read(io["coords"]).tables())
            else:
                table = next(snda.read(io["coords"]).tables())
                for column_name in table.column_names:
                    if column_name in io["dims"]:
                        assert table[column_name].to_pylist() == io["dims"][column_name]

        else:
            if io["throws"] is not None:
                with pytest.raises(io["throws"]):
                    r = snda.read(io["coords"])
                    if read_format == "coo":
                        next(r.coos())
                    elif read_format == "table":
                        next(r.tables())
            else:
                r = snda.read(io["coords"])
                if read_format == "coo":
                    tensor = next(r.coos())
                elif read_format == "table":
                    tensor = next(r.tables())
                assert tensor.shape == io["shape"]

        bad = False
        try:
            # attempt to write snda opened in read-only mode should fail
            snda.write(arrow_tensor)
            bad = True
        except Exception:
            pass
        assert not bad


@pytest.mark.parametrize(
    ("result_order", "want"),
    [
        (
            soma.ResultOrder.ROW_MAJOR,
            {
                "soma_dim_0": [2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4],
                "soma_dim_1": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5],
            },
        ),
        (
            "column-major",
            {
                "soma_dim_0": [2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4],
                "soma_dim_1": [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5],
            },
        ),
    ],
)
def test_result_order(
    tmp_path: pathlib.Path, result_order, want: Dict[str, List[float]]
):
    arrow_tensor = create_random_tensor("table", (5, 7), np.float32(), density=1)

    with soma.SparseNDArray.create(
        tmp_path.as_uri(), type=pa.float64(), shape=(5, 7)
    ) as write_arr:
        write_arr.write(arrow_tensor)
    with soma.open(tmp_path.as_uri()) as read_arr:
        assert isinstance(read_arr, soma.SparseNDArray)
        table = next(
            read_arr.read(
                [slice(2, 4), slice(1, 5)], result_order=result_order
            ).tables()
        )
        for col, values in want.items():
            assert table[col].to_pylist() == values


def test_sparse_nd_array_not_implemented(tmp_path):
    """Poke all of the expected not implemented API"""
    soma.SparseNDArray.create(
        tmp_path.as_posix(), type=pa.uint32(), shape=(99,)
    ).close()

    with soma.SparseNDArray.open(tmp_path.as_posix()) as a:
        with pytest.raises(NotImplementedError):
            next(a.read().dense_tensors())


def test_sparse_nd_array_error_corners(tmp_path):
    """Poke edge error handling"""
    with soma.SparseNDArray.create(
        tmp_path.as_posix(), type=pa.uint32(), shape=(99,)
    ) as a:
        a.write(
            create_random_tensor(
                format="coo", shape=(99,), dtype=np.uint32, density=0.1
            )
        )

        # Write should reject unknown types
        with pytest.raises(TypeError):
            a.write(pa.array(np.zeros((99,), dtype=np.uint32)))
        with pytest.raises(TypeError):
            a.write(pa.chunked_array([np.zeros((99,), dtype=np.uint32)]))

        # Write should reject wrong dimensionality
        with pytest.raises(ValueError):
            a.write(
                pa.SparseCSRMatrix.from_scipy(
                    sparse.random(10, 10, format="csr", dtype=np.uint32)
                )
            )
        with pytest.raises(ValueError):
            a.write(
                pa.SparseCSCMatrix.from_scipy(
                    sparse.random(10, 10, format="csc", dtype=np.uint32)
                )
            )

    with soma.SparseNDArray.open(tmp_path.as_posix()) as a:
        # other coord types are illegal
        with pytest.raises(TypeError):
            next(a.read("hi").tables())


@pytest.mark.parametrize(
    "bad_coords",
    [
        (slice(1, 10, 1),),  # explicit step
        (slice(32, 1),),  # stop < start
        (slice(-10, -2),),  # negative start & stop
        (slice(-32),),  # stop < entire domain
        (slice(10, -2),),  # stop < start
        (slice(100, None),),  # entire domain < start
        (slice(150, 200),),  # entire domain < start
        (slice(None), slice(None)),  # too many dims
    ],
)
def test_bad_coords(tmp_path, bad_coords):
    """
    Most illegal coords raise ValueError -- test for those.
    Oddly, some raise TypeError, which is covered in another
    test.
    """

    with soma.SparseNDArray.create(
        uri=tmp_path.as_posix(), type=pa.uint32(), shape=(99,)
    ) as a:
        a.write(
            create_random_tensor(
                format="coo", shape=(99,), dtype=np.uint32, density=0.1
            )
        )

    with _factory.open(tmp_path.as_posix()) as a:
        with pytest.raises(ValueError):
            next(a.read(bad_coords).tables())


def test_tile_extents(tmp_path):
    soma.SparseNDArray.create(
        tmp_path.as_posix(),
        type=pa.float32(),
        shape=(100, 10000),
        platform_config={
            "tiledb": {
                "create": {
                    "dims": {
                        "soma_dim_0": {"tile": 2048},
                        "soma_dim_1": {"tile": 2048},
                    }
                }
            }
        },
    ).close()

    with tiledb.open(tmp_path.as_posix()) as A:
        assert A.schema.domain.dim(0).tile == 100
        assert A.schema.domain.dim(1).tile == 2048


@pytest.mark.parametrize(
    "create_options,expected_schema_fields",
    (
        (
            {"allows_duplicates": True},
            {
                "validity_filters": tiledb.FilterList([tiledb.RleFilter()]),
                "allows_duplicates": True,
            },
        ),
        (
            {"allows_duplicates": False},
            {
                "validity_filters": tiledb.FilterList([tiledb.RleFilter()]),
                "allows_duplicates": False,
            },
        ),
        (
            {"validity_filters": ["NoOpFilter"], "allows_duplicates": False},
            {
                "validity_filters": tiledb.FilterList([tiledb.NoOpFilter()]),
                "allows_duplicates": False,
            },
        ),
    ),
)
def test_create_platform_config_overrides(
    tmp_path, create_options, expected_schema_fields
):
    uri = tmp_path.as_posix()
    soma.SparseNDArray.create(
        uri,
        type=pa.float64(),
        shape=(100, 100),
        platform_config={"tiledb": {"create": {**create_options}}},
    ).close()
    with tiledb.open(uri) as D:
        for k, v in expected_schema_fields.items():
            assert getattr(D.schema, k) == v


def test_timestamped_ops(tmp_path):
    # 2x2 array
    with soma.SparseNDArray.create(
        tmp_path.as_posix(),
        type=pa.uint16(),
        shape=(2, 2),
        tiledb_timestamp=10,
    ) as a:
        # write 1 into top-left entry @ t=10
        a.write(
            pa.SparseCOOTensor.from_scipy(
                sparse.coo_matrix(([1], ([0], [0])), shape=a.shape)
            )
        )

    # write 1 into bottom-right entry @ t=20
    with soma.SparseNDArray.open(
        tmp_path.as_posix(), mode="w", tiledb_timestamp=20
    ) as a:
        a.write(
            pa.SparseCOOTensor.from_scipy(
                sparse.coo_matrix(([1], ([1], [1])), shape=a.shape)
            )
        )

    # read with no timestamp args & see both 1s
    with soma.SparseNDArray.open(tmp_path.as_posix()) as a:
        assert a.read().coos().concat().to_scipy().todense().tolist() == [
            [1, 0],
            [0, 1],
        ]
        assert a.nnz == 2

    # read @ t=15 & see only the first write
    with soma.SparseNDArray.open(
        tmp_path.as_posix(), context=SOMATileDBContext(timestamp=15)
    ) as a:
        assert a.read().coos().concat().to_scipy().todense().tolist() == [
            [1, 0],
            [0, 0],
        ]
        assert a.nnz == 1


def test_empty_indexed_read(tmp_path):
    """
    Verify that queries expected to return empty results actually
    work. There are edge cases around SparseTensors, which are unable
    to represent empty arrays.
    """
    shape = (10, 100)
    soma.SparseNDArray.create(
        tmp_path.as_posix(), type=pa.uint16(), shape=shape
    ).close()

    data = create_random_tensor("coo", shape, np.float64, 1.0)
    with soma.SparseNDArray.open(tmp_path.as_posix(), "w") as a:
        a.write(data)

    with soma.SparseNDArray.open(tmp_path.as_posix()) as a:
        coords = [slice(None), slice(None)]
        assert sum(len(t) for t in a.read(coords).tables()) == 1000

        coords = [[3], [4]]
        assert sum(len(t) for t in a.read(coords).tables()) == 1

        coords = [[3], []]
        assert sum(len(t) for t in a.read(coords).tables()) == 0

        coords = [[], [4]]
        assert sum(len(t) for t in a.read(coords).tables()) == 0


@pytest.fixture
def a_random_sparse_nd_array(tmp_path, shape: Tuple[int, ...], density: float) -> str:
    uri = tmp_path.as_posix()
    dtype = np.float32
    with soma.SparseNDArray.create(
        uri, type=pa.from_numpy_dtype(dtype), shape=shape
    ) as a:
        a.write(create_random_tensor("table", shape, dtype, density))
    return uri


@pytest.mark.parametrize(
    # specify coords using SOMA semantics, ie. closed range
    # use density to keep the sparse array size "reasonable"
    "density,shape,coords",
    [
        (0.1, (1_000, 100), ()),
        (0.1, (1_000, 100), (None,)),
        (0.1, (1_000, 100), (slice(None),)),
        (0.1, (1_000, 100), (slice(10, 20),)),
        (0.1, (1_000, 100), (None, slice(10, 20))),
        (0.1, (1_000, 100), (slice(10, 20), slice(2, 33))),
        (0.1, (1_000, 100), (1,)),
        (0.1, (1_000, 100), (None, 3)),
        (0.1, (1_000, 100), (10, 3)),
        (0.1, (1_000, 100), (slice(None), 1)),
        (0.1, (1_000, 100), ([3, 99, 101, 102, 77, 0],)),
        (0.1, (1_000, 100), (slice(500), [3, 99, 33, 77, 0])),
        (0.1, (1_000, 100), (np.arange(0, 101),)),
        (0.1, (1_000, 100), (np.arange(0, 1000), [0, 10, 20])),
        (0.1, (1_000, 100), (np.arange(0, 1000), slice(2, 99))),
        (0.01, (10_101, 389), ()),
        (0.01, (10_101, 389), (None,)),
        (0.01, (10_101, 389), (slice(None),)),
        (0.01, (10_101, 389), (slice(33, 1048),)),
        (0.01, (10_101, 389), (slice(33, 1048), slice(45, 333))),
    ],
)
@pytest.mark.parametrize("step", [9, 77, 243, 1001, 2**16])
def test_scipy_iter_reindexed(
    a_random_sparse_nd_array: str, coords: Tuple[Any, ...], step: int
) -> None:
    """
    Verify that simple use of scipy iterator works.
    """

    def _slice_sp(
        coo: sparse.coo_matrix, _coords: Tuple[Any, ...]
    ) -> sparse.coo_matrix:
        """
        Slice from the COO, accomodating conversion from closed range to half-open range
        slices, plus quirks of scipy.sparse which can't slice on multiple dimensions in
        all cases (but handles it fine if you slice one dim at a time).
        """
        csr = coo.tocsr()
        for i, c in enumerate(_coords):
            if isinstance(c, slice) and c.stop is not None:
                c = slice(c.start, c.stop + 1)
            if c is None:
                c = slice(None)
            _coord = tuple([slice(None)] * (i) + [c])
            csr = operator.getitem(csr, _coord)
        return csr.tocoo()

    # these are not pytest params to speed up tests (by reducing the number of SOMA arrays created)
    for axis, compress, reindex_sparse_axis in [
        (a, c, ri) for a in (1, 0) for c in (False, True) for ri in (True, False)
    ]:
        minor_axis = 1 - axis
        with soma.open(a_random_sparse_nd_array, mode="r") as A:
            truth_coo_tbl = A.read().tables().concat()
            truth_coo = sparse.coo_matrix(
                (
                    truth_coo_tbl.column(2).to_numpy(),
                    (
                        truth_coo_tbl.column(0).to_numpy(),
                        truth_coo_tbl.column(1).to_numpy(),
                    ),
                ),
                shape=A.shape,
            )
            truth_coo = _slice_sp(truth_coo, coords)

            results = []
            for joinids, sp in A.read(coords).scipy(
                axis=axis,
                compress=compress,
                step=step,
                reindex_sparse_axis=reindex_sparse_axis,
            ):
                # check for expected type
                assert len(joinids) == 2 and isinstance(joinids, tuple)
                assert all(isinstance(j, np.ndarray) for j in joinids)
                if not compress:
                    assert isinstance(sp, sparse.coo_matrix)
                elif axis == 0:
                    assert isinstance(sp, sparse.csr_matrix)
                else:
                    assert isinstance(sp, sparse.csc_matrix)

                # check for expected shape
                assert len(joinids[axis]) == min(step, sp.shape[axis])
                assert (
                    not reindex_sparse_axis
                    or len(joinids[minor_axis]) == sp.shape[minor_axis]
                )

                # internal layout (dups, ordering, etc)
                assert sp.has_canonical_format
                if compress:
                    sp.check_format(full_check=True)

                # sanity check coordinates
                if not reindex_sparse_axis:
                    if not compress:
                        minor = sp.col if axis == 0 else sp.row
                    else:
                        minor = sp.indices
                    assert np.isin(minor, joinids[minor_axis]).all()

                results.append(sp)

            # check vs ground truth - only implemented if reindex_sparse_axis == True
            if reindex_sparse_axis:
                stacked = (
                    sparse.vstack(results) if axis == 0 else sparse.hstack(results)
                )
                assert (
                    truth_coo.dtype == stacked.dtype
                    and truth_coo.shape == stacked.shape
                    and (truth_coo != stacked.tocoo()).nnz == 0
                )


@pytest.mark.parametrize("density,shape", [(0.1, (100, 100))])
def test_scipy_iter_error_checks(
    a_random_sparse_nd_array: str, shape: Tuple[int, ...]
) -> None:
    with soma.open(a_random_sparse_nd_array, mode="r") as A:
        with pytest.raises(soma.SOMAError):
            A.read().scipy(axis=2)

        with pytest.raises(ValueError):
            A.read(result_order=soma.ResultOrder.ROW_MAJOR).scipy(axis=1)

        with pytest.raises(ValueError):
            A.read(result_order=soma.ResultOrder.COLUMN_MAJOR).scipy(axis=0)


@pytest.mark.parametrize("density,shape", [(0.1, (4, 8, 16))])
def test_scipy_iter_not_2D(
    a_random_sparse_nd_array: str, shape: Tuple[int, ...]
) -> None:
    with soma.open(a_random_sparse_nd_array, mode="r") as A:
        with pytest.raises(soma.SOMAError):
            A.read().scipy()
