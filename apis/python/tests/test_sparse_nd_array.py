from typing import Optional, Tuple, Union

import numpy as np
import pyarrow as pa
import pytest
import scipy.sparse as sparse

import tiledbsoma as soma

from . import NDARRAY_ARROW_TYPES_NOT_SUPPORTED, NDARRAY_ARROW_TYPES_SUPPORTED

AnySparseTensor = Union[pa.SparseCOOTensor, pa.SparseCSRMatrix, pa.SparseCSCMatrix]


def test_sparse_nd_array_ok_no_storage():
    arr = soma.SparseNdArray(uri="/foo/bar/")
    assert arr.uri == "/foo/bar/"
    assert not arr.exists()
    assert arr.soma_type == "SOMASparseNdArray"


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

    a = soma.SparseNdArray(uri=tmp_path.as_posix())
    a.create(element_type, shape)
    assert a.soma_type == "SOMASparseNdArray"
    assert a.uri == tmp_path.as_posix()
    assert a.ndim == len(shape)
    assert a.shape == tuple(shape)
    assert a.is_sparse is True
    assert a.exists()

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
    a = soma.SparseNdArray(uri=tmp_path.as_posix())
    with pytest.raises(TypeError):
        a.create(element_type, shape)
    assert not a.exists()


def test_sparse_nd_array_delete(tmp_path):
    a = soma.SparseNdArray(uri=tmp_path.as_posix())
    a.create(pa.int8(), (100, 100))
    assert a.exists()

    a.delete()
    assert not a.exists()
    assert not tmp_path.exists()  # check underlying file system

    # should be silent about non-existent object
    assert a.delete() is None
    assert soma.SparseNdArray(uri="no such array").delete() is None


def create_random_tensor(
    format: str,
    shape: Tuple[int, ...],
    dtype: np.dtype,
    density: Optional[float] = 0.33,
):
    """
    Create a random tensor/table of specified format, shape and dtype.

    Guarantee: there will be NO duplicate values in the tensor, which makes
    it simpler to validate (see `tensors_are_same_value`)
    """
    rng = np.random.default_rng()
    ndim = len(shape)
    assert density > 0 and density <= 1

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
        ((10, 21), "csr"),
        ((10, 21), "csc"),
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
    assert format in ("coo", "csr", "csc")
    assert not (format in ("csc", "csr") and len(shape) != 2)

    a = soma.SparseNdArray(tmp_path.as_posix())
    a.create(pa.float64(), shape)
    assert a.exists()
    assert a.shape == shape

    # Make a random sample in the desired format
    # As discussed in the SparseNdArray implementation, Arrow SparseTensor objects can't be zero-length
    # so we must be prepared for StopIteration on reading them. It simplifies unit-test logic to use
    # occupation density of 1.0 for this test.
    data = create_random_tensor(format, shape, np.float64, 1.0)
    a.write_sparse_tensor(data)
    del a

    # Read back and validate
    b = soma.SparseNdArray(tmp_path.as_posix())

    t = next(b.read_sparse_tensor((slice(None),) * len(shape), format=format))
    assert tensors_are_same_value(t, data)

    t = next(b.read_sparse_tensor((0,) * len(shape), format=format))
    assert t.shape == shape


@pytest.mark.parametrize("shape", [(10,), (23, 4), (5, 3, 1), (8, 4, 2, 30)])
@pytest.mark.parametrize("test_enumeration", range(10))
def test_sparse_nd_array_read_write_table(
    tmp_path, shape: Tuple[int, ...], test_enumeration: int
):
    a = soma.SparseNdArray(tmp_path.as_posix())
    a.create(pa.float32(), shape)
    assert a.exists()
    assert a.shape == shape

    # make a random sample in the desired format
    data = create_random_tensor("table", shape, np.float32)
    a.write_table(data)
    del a

    # Read back and validate
    b = soma.SparseNdArray(tmp_path.as_posix())
    t = next(b.read_table((slice(None),) * len(shape)))
    assert isinstance(t, pa.Table)
    assert tables_are_same_value(data, t)


@pytest.mark.parametrize("dtype", [np.float32, np.float64, np.int32, np.int64])
@pytest.mark.parametrize("shape", [(1,), (23, 14), (35, 3, 2), (8, 4, 2, 30)])
def test_sparse_nd_array_read_as_pandas(
    tmp_path, dtype: np.dtype, shape: Tuple[int, ...]
):

    dtype = np.dtype(dtype)
    a = soma.SparseNdArray(tmp_path.as_posix())
    a.create(pa.from_numpy_dtype(dtype), shape)
    assert a.exists()
    assert a.shape == shape

    # make a random sample in the desired format
    data = create_random_tensor("table", shape, dtype)
    a.write_table(data)

    df = a.read_as_pandas_all()

    dim_names = [f"soma_dim_{n}" for n in range(len(shape))]
    assert df.sort_values(by=dim_names, ignore_index=True).equals(
        data.to_pandas().sort_values(by=dim_names, ignore_index=True)
    )


def test_empty_read(tmp_path):
    """
    Verify that queries expected to return empty results actually
    work. There are edge cases around SparseTensors, which are unable
    to represent empty arrays.
    """
    a = soma.SparseNdArray(uri=tmp_path.as_posix())
    a.create(type=pa.uint16(), shape=(10, 100))
    assert a.exists()

    #
    # First, test reads of zero element array
    #

    # These work as expected
    coords = (slice(None),)
    assert sum(len(t) for t in a.read_table(coords)) == 0
    assert (
        sum(t.non_zero_length for t in a.read_sparse_tensor(coords, format="coo")) == 0
    )
    assert (
        sum(t.non_zero_length for t in a.read_sparse_tensor(coords, format="csr")) == 0
    )
    assert (
        sum(t.non_zero_length for t in a.read_sparse_tensor(coords, format="csc")) == 0
    )

    #
    # Next, test empty queries on non-empty array
    #
    a.write_sparse_tensor(
        pa.SparseCOOTensor.from_scipy(
            sparse.coo_matrix(([1], ([0], [0])), shape=a.shape)
        )
    )
    assert sum(len(t) for t in a.read_table((slice(None),))) == 1

    coords = (1, 1)  # no element at this coordinate
    assert sum(len(t) for t in a.read_table(coords)) == 0
    assert (
        sum(t.non_zero_length for t in a.read_sparse_tensor(coords, format="coo")) == 0
    )
    assert (
        sum(t.non_zero_length for t in a.read_sparse_tensor(coords, format="csr")) == 0
    )
    assert (
        sum(t.non_zero_length for t in a.read_sparse_tensor(coords, format="csc")) == 0
    )


@pytest.mark.parametrize("shape", [(), (0,), (10, 0), (0, 10), (1, 2, 0)])
def test_zero_length_fail(tmp_path, shape):
    """Zero length dimensions are expected to fail"""
    a = soma.SparseNdArray(tmp_path.as_posix())
    with pytest.raises(ValueError):
        a.create(type=pa.float32(), shape=shape)


def test_sparse_nd_array_nnz(tmp_path):
    """
    This operation is currently unimplemented. This test is
    designed to start failing as soon as it is implemented,
    to provide a reminder to create a real test.

    Just remove the pytest.raises when implemented
    """
    a = soma.SparseNdArray(tmp_path.as_posix())
    a.create(type=pa.int32(), shape=(10, 10, 10))
    assert a.nnz == 0

    t: pa.SparseCOOTensor = create_random_tensor(
        "coo", a.shape, pa.int32().to_pandas_dtype(), 0.1
    )
    a.write_sparse_tensor(t)
    assert t.non_zero_length == a.nnz


def test_sparse_nd_array_reshape(tmp_path):
    """
    Reshape currently unimplemented.
    """
    a = soma.SparseNdArray(tmp_path.as_posix())
    a.create(type=pa.int32(), shape=(10, 10, 10))
    with pytest.raises(NotImplementedError):
        assert a.reshape((100, 10, 1))


@pytest.mark.parametrize(
    "io",
    [

        # Coords is None
        {
            "shape": (4,),
            "coords": None,
            "dims": {
                "soma_dim_0": [0,1,2,3],
            },
            "throws": None,
        },

        # Coords has None in a slot
        {
            "shape": (4,),
            "coords": (None,),
            "dims": {
                "soma_dim_0": [0,1,2,3],
            },
            "throws": None,
        },

        # Coords has int in a slot
        {
            "shape": (4,),
            "coords": (1,),
            "dims": {
                "soma_dim_0": [1],
            },
            "throws": None,
        },
        {
            "shape": (4, 6),
            "coords": (0, 0),
            "dims": {
                "soma_dim_0": [0],
                "soma_dim_1": [0],
            },
            "throws": None,
        },

        # Coords doesn't specify all dimensions, so the rest are implicit-all
        {
            "shape": (4, 6),
            "coords": (0,),
            "dims": {
                "soma_dim_0": [0, 0, 0, 0, 0, 0],
                "soma_dim_1": [0, 1, 2, 3, 4, 5],
            },
            "throws": None,
        },

        # Coords specifies too many dimensions
        {
            "shape": (4, 6),
            "coords": (0,0,0),
            "dims": {
                "soma_dim_0": [0, 0, 0, 0, 0, 0],
                "soma_dim_1": [0, 1, 2, 3, 4, 5],
            },
            "throws": ValueError,
        },

        {
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
            "shape": (4, 6),
            "coords": (3, 4),
            "dims": {
                "soma_dim_0": [3],
                "soma_dim_1": [4],
            },
            "throws": None,
        },

        {
            "shape": (4, 6),
            "coords": (slice(1, 2), slice(3, 4)),
            "dims": {
                "soma_dim_0": [1, 1, 2, 2],
                "soma_dim_1": [3, 4, 3, 4],
            },
            "throws": None,
        },

        {
            "shape": (4, 6),
            "coords": (slice(None), slice(3, 4)),
            "dims": {
                "soma_dim_0": [0, 0, 1, 1, 2, 2, 3, 3],
                "soma_dim_1": [3, 4, 3, 4, 3, 4, 3, 4],
            },
            "throws": None,
        },

        {
            "shape": (4, 6),
            "coords": (slice(1, 2), slice(None)),
            "dims": {
                "soma_dim_0": [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2],
                "soma_dim_1": [0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5],
            },
            "throws": None,
        },

        {
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

    ],
)
def test_sparse_nd_array_table_slicing(tmp_path, io):

    # Set up contents
    pacoo = create_random_tensor(
        format="coo",
        shape=io["shape"],
        dtype=np.float32(),
        density=1.0,
    )

    snda = soma.SparseNdArray(tmp_path.as_posix())
    snda.create(pa.float64(), io["shape"])
    snda.write_sparse_tensor(pacoo)

    if io["throws"] is not None:
        with pytest.raises(io["throws"]):
            next(snda.read_table(io["coords"]))
    else:
        table = next(snda.read_table(io["coords"]))
        for column_name in table.column_names:
            if column_name in io["dims"]:
                assert table[column_name].to_pylist() == io["dims"][column_name]

    # TODO:
    # * read methods:
    #   o snda.read_as_pandas -> just a wrapper around read_table -> to_pandas
    #   o snda.read_sparse_tensor -> reshape around read_table
    # * ndim
    # * format
    #   o coo require 2d
    #   o nonesuch format -> NotImplementedError
    # * StopIteration on empty query
    # * ids type-switching borrow from 470
