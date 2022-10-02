from typing import Tuple, Union

import numpy as np
import pyarrow as pa
import pytest
import scipy.sparse as sparse

import tiledbsoma as soma

from . import NDARRAY_ARROW_TYPES_NOT_SUPPORTED, NDARRAY_ARROW_TYPES_SUPPORTED

"""
comment will be removed when test rework is complete
TODO:
- [X] create
- [X] delete
- [X] exists
- [X] get type
- [X] get shape
- [X] get ndims
- [X] get schema
- [X] get is_sparse
- [X] get metadata
- [ ] get nnz
- [X] read
- [X] write
- [ ] reshape
"""
AnySparseTensor = Union[pa.SparseCOOTensor, pa.SparseCSRMatrix, pa.SparseCSCMatrix]


def test_soma_sparse_nd_array_ok_no_storage():
    arr = soma.SOMASparseNdArray(uri="/foo/bar")
    assert arr.uri == "/foo/bar"
    assert arr.name == "bar"
    assert not arr.exists()
    assert arr.type == "SOMASparseNdArray"


@pytest.mark.parametrize(
    "shape", [(10,), (1, 100), (10, 1, 100), (2, 4, 6, 8), [1], (1, 2, 3, 4, 5)]
)
@pytest.mark.parametrize("element_type", NDARRAY_ARROW_TYPES_SUPPORTED)
def test_soma_sparse_nd_array_create_ok(
    tmp_path, shape: Tuple[int, ...], element_type: pa.DataType
):
    """
    Test all cases we expect "create" to succeed.
    """
    assert pa.types.is_primitive(element_type)  # sanity check incoming params

    a = soma.SOMASparseNdArray(uri=tmp_path.as_posix())
    a.create(element_type, shape)
    assert a.type == "SOMASparseNdArray"
    assert a.uri == tmp_path.as_posix()
    assert a.ndims == len(shape)
    assert a.shape == tuple(shape)
    assert a.is_sparse is True
    assert a.exists()

    assert a.schema is not None
    expected_field_names = ["data"] + [f"__dim_{d}" for d in range(len(shape))]
    assert set(a.schema.names) == set(expected_field_names)
    for d in range(len(shape)):
        assert a.schema.field(f"__dim_{d}").type == pa.uint64()
    assert a.schema.field("data").type == element_type


@pytest.mark.parametrize("shape", [(10,)])
@pytest.mark.parametrize("element_type", NDARRAY_ARROW_TYPES_NOT_SUPPORTED)
def test_soma_sparse_nd_array_create_fail(
    tmp_path, shape: Tuple[int, ...], element_type: pa.DataType
):
    a = soma.SOMASparseNdArray(uri=tmp_path.as_posix())
    with pytest.raises(TypeError):
        a.create(element_type, shape)
    assert not a.exists()


def test_soma_sparse_nd_array_delete(tmp_path):
    a = soma.SOMASparseNdArray(uri=tmp_path.as_posix())
    a.create(pa.int8(), (100, 100))
    assert a.exists()

    a.delete()
    assert not a.exists()
    assert not tmp_path.exists()  # check underlying file system

    # should be silent about non-existent object
    assert a.delete() is None
    assert soma.SOMASparseNdArray(uri="no such array").delete() is None


def create_random_tensor(format: str, shape: Tuple[int, ...], dtype: np.dtype):
    """
    Create a random tensor/table of specified format, shape and dtype.

    Guarantee: there will be NO duplicate values in the tensor, which makes
    it simpler to validate (see `tensors_are_same_value`)
    """
    rng = np.random.default_rng()
    ndim = len(shape)

    assert format in ["coo", "csc", "csr", "table"], "Unimplemented format"

    if format == "coo":
        nrec = rng.integers(np.prod(shape)) + 1
        data = rng.choice(10 * nrec, size=nrec, replace=False).astype(dtype)
        all_coords = np.array(
            np.meshgrid(*tuple(np.arange(dim_len) for dim_len in shape))
        ).T.reshape(-1, ndim)
        coords = rng.choice(all_coords, nrec, replace=False)
        return pa.SparseCOOTensor.from_numpy(data, coords, shape=shape)

    if format == "table":
        nrec = rng.integers(np.prod(shape)) + 1
        data = rng.choice(10 * nrec, size=nrec, replace=False).astype(dtype)
        all_coords = np.array(
            np.meshgrid(*tuple(np.arange(dim_len) for dim_len in shape))
        ).T.reshape(-1, ndim)
        coords = rng.choice(all_coords, nrec, replace=False).T
        pydict = {
            f"__dim_{n}": pa.array(coords[n], type=pa.uint64()) for n in range(ndim)
        }
        pydict.update({"data": pa.array(data)})
        return pa.Table.from_pydict(pydict)

    if format == "csc":
        assert ndim == 2
        return pa.SparseCSCMatrix.from_scipy(
            sparse.random(
                *shape,
                density=0.33,
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
                density=0.33,
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
    if a.field("data").type != b.field("data").type:
        return False
    for tbl in (a, b):
        if not all(
            tbl.field(f"__dim_{n}").type == pa.uint64()
            for n in range(tbl.num_columns - 1)
        ):
            return False

    ndim = a.shape[1] - 1
    ai = a.column("data").to_numpy().argsort()
    bi = b.column("data").to_numpy().argsort()
    if not np.array_equal(
        np.take_along_axis(a.column("data").to_numpy(), ai, axis=0),
        np.take_along_axis(b.column("data").to_numpy(), bi, axis=0),
    ):
        return False

    for n in range(ndim):
        dim_name = f"__dim_{n}"
        if (
            a.field(dim_name).type != pa.uint64()
            or b.field(dim_name).type != pa.uint64()
        ):
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
        ((10, 21), "csr"),
        ((10, 21), "csc"),
        ((10, 20, 2), "coo"),
        ((2, 4, 6, 8), "coo"),
        ((1, 2, 4, 6, 8), "coo"),
    ],
)
@pytest.mark.parametrize("test_enumeration", range(10))
def test_soma_sparse_nd_array_read_write_sparse_tensor(
    tmp_path,
    shape: Tuple[int, ...],
    format: str,
    test_enumeration: int,
):
    # Test sanity: Tensor only, and CSC and CSR only support 2D, so fail any nonsense configs
    assert format in ("coo", "csr", "csc")
    assert not (format in ("csc", "csr") and len(shape) != 2)

    a = soma.SOMASparseNdArray(tmp_path.as_posix())
    a.create(pa.float64(), shape)
    assert a.exists()
    assert a.shape == shape

    # make a random sample in the desired format
    data = create_random_tensor(format, shape, np.float64)
    a.write_sparse_tensor(data)
    del a

    # Read back and validate
    b = soma.SOMASparseNdArray(tmp_path.as_posix())
    t = next(b.read_sparse_tensor((slice(None),) * len(shape), format=format))
    assert tensors_are_same_value(t, data)


@pytest.mark.parametrize("shape", [(10,), (23, 4), (5, 3, 1), (8, 4, 2, 30)])
@pytest.mark.parametrize("test_enumeration", range(10))
def test_soma_sparse_nd_array_read_write_table(
    tmp_path, shape: Tuple[int, ...], test_enumeration: int
):
    a = soma.SOMASparseNdArray(tmp_path.as_posix())
    a.create(pa.float32(), shape)
    assert a.exists()
    assert a.shape == shape

    # make a random sample in the desired format
    data = create_random_tensor("table", shape, np.float32)
    a.write_table(data)
    del a

    # Read back and validate
    b = soma.SOMASparseNdArray(tmp_path.as_posix())
    t = next(b.read_table((slice(None),) * len(shape)))
    assert isinstance(t, pa.Table)
    assert tables_are_same_value(data, t)


@pytest.mark.parametrize("dtype", [np.float32, np.float64, np.int32, np.int64])
@pytest.mark.parametrize("shape", [(1,), (23, 14), (35, 3, 2), (8, 4, 2, 30)])
def test_soma_sparse_nd_array_read_as_pandas(
    tmp_path, dtype: np.dtype, shape: Tuple[int, ...]
):

    dtype = np.dtype(dtype)
    a = soma.SOMASparseNdArray(tmp_path.as_posix())
    a.create(pa.from_numpy_dtype(dtype), shape)
    assert a.exists()
    assert a.shape == shape

    # make a random sample in the desired format
    data = create_random_tensor("table", shape, dtype)
    a.write_table(data)

    df = a.read_as_pandas_all()

    dim_names = [f"__dim_{n}" for n in range(len(shape))]
    assert df.sort_values(by=dim_names, ignore_index=True).equals(
        data.to_pandas().sort_values(by=dim_names, ignore_index=True)
    )


@pytest.mark.parametrize("shape", [(), (0,), (10, 0), (0, 10), (1, 2, 0)])
def test_zero_length_fail(tmp_path, shape):
    """Zero length dimensions are expected to fail"""
    a = soma.SOMASparseNdArray(tmp_path.as_posix())
    with pytest.raises(ValueError):
        a.create(type=pa.float32(), shape=shape)
