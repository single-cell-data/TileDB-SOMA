from __future__ import annotations

import contextlib
import datetime
import gc
import itertools
import json
import operator
import pathlib
import sys
from concurrent import futures
from typing import Any, Dict, List, Tuple, Union
from unittest import mock

import numpy as np
import pyarrow as pa
import pytest
import scipy.sparse as sparse

import tiledbsoma as soma
from tiledbsoma import _factory
from tiledbsoma.options import SOMATileDBContext

from . import NDARRAY_ARROW_TYPES_NOT_SUPPORTED, NDARRAY_ARROW_TYPES_SUPPORTED
from ._util import raises_no_typeguard

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

    with raises_no_typeguard(TypeError):
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

    # TODO: more testing with current-domain feature integrated
    # https://github.com/single-cell-data/TileDB-SOMA/issues/2407
    assert isinstance(a.maxshape, tuple)
    assert len(a.maxshape) == len(a.shape)
    for ms, s in zip(a.maxshape, a.shape):
        assert ms >= s

    assert a.is_sparse is True

    assert a.schema is not None
    expected_field_names = ["soma_data"] + [f"soma_dim_{d}" for d in range(len(shape))]
    assert set(a.schema.names) == set(expected_field_names)
    for d in range(len(shape)):
        assert a.schema.field(f"soma_dim_{d}").type == pa.int64()
    assert a.schema.field("soma_data").type == element_type
    assert not a.schema.field("soma_data").nullable

    # Check with open binding
    with contextlib.closing(
        soma.pytiledbsoma.SOMASparseNDArray.open(
            tmp_path.as_posix(),
            soma.pytiledbsoma.OpenMode.read,
            soma.pytiledbsoma.SOMAContext(),
        )
    ) as b:
        assert a.schema == b.schema

    # Ensure read mode uses clib object
    with soma.SparseNDArray.open(tmp_path.as_posix(), "r") as A:
        assert isinstance(A._handle._handle, soma.pytiledbsoma.SOMASparseNDArray)

    # Ensure write mode uses clib object
    with soma.SparseNDArray.open(tmp_path.as_posix(), "w") as A:
        assert isinstance(A._handle._handle, soma.pytiledbsoma.SOMASparseNDArray)

    # Ensure it cannot be opened by another type
    with pytest.raises(soma.SOMAError):
        soma.DataFrame.open(tmp_path.as_posix())

    with pytest.raises(soma.SOMAError):
        soma.DenseNDArray.open(tmp_path.as_posix())

    with pytest.raises(soma.SOMAError):
        soma.PointCloudDataFrame.open(tmp_path.as_posix())

    with pytest.raises(soma.SOMAError):
        soma.Collection.open(tmp_path.as_posix())

    with pytest.raises(soma.SOMAError):
        soma.Experiment.open(tmp_path.as_posix())

    with pytest.raises(soma.SOMAError):
        soma.Measurement.open(tmp_path.as_posix())

    with pytest.raises(soma.SOMAError):
        soma.Scene.open(tmp_path.as_posix())

    with pytest.raises(soma.SOMAError):
        soma.MultiscaleImage.open(tmp_path.as_posix())


def test_sparse_nd_array_reopen(tmp_path):
    soma.SparseNDArray.create(
        tmp_path.as_posix(), type=pa.float64(), shape=(1,), tiledb_timestamp=1
    )

    with soma.SparseNDArray.open(tmp_path.as_posix(), "r", tiledb_timestamp=1) as A1:
        with raises_no_typeguard(ValueError):
            A1.reopen("invalid")

        with A1.reopen("w", tiledb_timestamp=2) as A2:
            with A2.reopen("r", tiledb_timestamp=3) as A3:
                assert A1.mode == "r"
                assert A2.mode == "w"
                assert A3.mode == "r"
                assert A1.tiledb_timestamp_ms == 1
                assert A2.tiledb_timestamp_ms == 2
                assert A3.tiledb_timestamp_ms == 3

    ts1 = datetime.datetime(2023, 1, 1, 1, 0, tzinfo=datetime.timezone.utc)
    ts2 = datetime.datetime(2024, 1, 1, 1, 0, tzinfo=datetime.timezone.utc)
    with soma.SparseNDArray.open(tmp_path.as_posix(), "r", tiledb_timestamp=ts1) as A1:
        with A1.reopen("r", tiledb_timestamp=ts2) as A2:
            assert A1.mode == "r"
            assert A2.mode == "r"
            assert A1.tiledb_timestamp == ts1
            assert A2.tiledb_timestamp == ts2

    with soma.SparseNDArray.open(tmp_path.as_posix(), "w") as A1:
        with A1.reopen("w", tiledb_timestamp=None) as A2:
            with A2.reopen("w") as A3:
                assert A1.mode == "w"
                assert A2.mode == "w"
                assert A3.mode == "w"
                now = datetime.datetime.now(datetime.timezone.utc)
                assert A1.tiledb_timestamp <= now
                assert A2.tiledb_timestamp <= now
                assert A3.tiledb_timestamp <= now


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
    if type(a) is not type(b):
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
    uri = tmp_path.as_posix()

    # Test sanity: Tensor only, and CSC and CSR only support 2D, so fail any nonsense configs
    assert format == "coo"

    a = soma.SparseNDArray.create(uri, type=pa.float64(), shape=shape)
    assert a.shape == shape

    # Make a random sample in the desired format
    # As discussed in the SparseNDArray implementation, Arrow SparseTensor objects can't be zero-length
    # so we must be prepared for StopIteration on reading them. It simplifies unit-test logic to use
    # occupation density of 1.0 for this test.
    data = create_random_tensor(format, shape, np.float64, 1.0)
    with raises_no_typeguard(TypeError):
        # non-arrow write
        a.write(data.to_numpy())
    a.write(data)
    a.close()
    del a

    # Array write should fail if array opened in read mode
    with soma.SparseNDArray.open(uri) as a:
        with pytest.raises(soma.SOMAError):
            a.write(data)

    # Read back and validate
    with soma.SparseNDArray.open(uri) as b:
        t = b.read((slice(None),) * len(shape)).coos().concat()

        assert tensors_are_same_value(t, data)

        t = next(b.read((0,) * len(shape)).coos())

        assert t.shape == shape

    with soma.SparseNDArray.open(uri) as A:
        assert A.is_sparse
        assert not A.schema_config_options().allows_duplicates


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

    with soma.SparseNDArray.open(tmp_path.as_posix()) as A:
        assert A.is_sparse
        assert not A.schema_config_options().allows_duplicates


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

    with soma.SparseNDArray.open(tmp_path.as_posix()) as A:
        assert A.is_sparse
        assert not A.schema_config_options().allows_duplicates


@pytest.mark.parametrize("shape_is_nones", [True, False])
@pytest.mark.parametrize("element_type", NDARRAY_ARROW_TYPES_SUPPORTED)
def test_sparse_nd_array_shaping(tmp_path, shape_is_nones, element_type):
    uri = tmp_path.as_posix()

    shape = [2, 3]

    soma.SparseNDArray.create(
        uri,
        type=element_type,
        shape=shape,
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
            snda.resize([3, 3])
        with soma.SparseNDArray.open(uri, "w") as snda:
            snda.write(batch2)
    else:
        # tiledbsoma._exception.SOMAError: [TileDB::Dimension] Error:
        # Coordinate 2 is out of domain bounds [0, 1] on dimension 'soma_dim_0'
        with pytest.raises(soma.SOMAError):
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
            "throws": (soma.SOMAError),
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
        {
            "name": "2D coords=(pa_chunked[1, 2], pa64[3, 4])",
            "shape": (9, 11),
            "coords": (
                pa.chunked_array([[1], [2]]),
                pa.array([3, 4]),
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
    uri = tmp_path.as_posix()

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

    with soma.SparseNDArray.create(uri, type=pa.float64(), shape=io["shape"]) as snda_w:
        snda_w.write(arrow_tensor)

    with soma.SparseNDArray.open(uri) as snda:
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
        with raises_no_typeguard(TypeError):
            a.write(pa.array(np.zeros((99,), dtype=np.uint32)))
        with raises_no_typeguard(TypeError):
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
        with raises_no_typeguard(TypeError):
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
    uri = tmp_path.as_posix()
    soma.SparseNDArray.create(
        uri,
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

    with soma.SparseNDArray.open(tmp_path.as_posix()) as A:
        dim_info = json.loads(A.schema_config_options().dims)
        assert int(dim_info["soma_dim_0"]["tile"]) == 2048
        assert int(dim_info["soma_dim_1"]["tile"]) == 2048


@pytest.mark.parametrize(
    "create_options,expected_schema_fields",
    (
        (
            {"allows_duplicates": True},
            {
                "validity_filters": [{"COMPRESSION_LEVEL": -1, "name": "RLE"}],
                "allows_duplicates": True,
            },
        ),
        (
            {"allows_duplicates": False},
            {
                "validity_filters": [{"COMPRESSION_LEVEL": -1, "name": "RLE"}],
                "allows_duplicates": False,
            },
        ),
        (
            {"validity_filters": ["NoOpFilter"], "allows_duplicates": False},
            {
                "validity_filters": [{"name": "NOOP"}],
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

    with soma.SparseNDArray.open(tmp_path.as_posix()) as A:
        cfg = A.schema_config_options()
        assert expected_schema_fields["validity_filters"] == json.loads(
            cfg.validity_filters
        )
        assert expected_schema_fields["allows_duplicates"] == cfg.allows_duplicates


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
def a_soma_context() -> SOMATileDBContext:
    return SOMATileDBContext(
        tiledb_config={
            "soma.init_buffer_bytes": 128 * 1024**2,
            "tiledb.init_buffer_bytes": 128 * 1024**2,
        }
    )


@pytest.fixture
def a_random_sparse_nd_array(
    tmp_path, a_soma_context: SOMATileDBContext, shape: Tuple[int, ...], density: float
) -> str:
    uri = tmp_path.as_posix()
    dtype = np.float32
    with soma.SparseNDArray.create(
        uri, type=pa.from_numpy_dtype(dtype), shape=shape, context=a_soma_context
    ) as a:
        a.write(create_random_tensor("table", shape, dtype, density))
    return uri


@pytest.mark.parametrize(
    # specify coords using SOMA semantics, ie. closed range
    # use density to keep the sparse array memory use "reasonable"
    "density,shape,coords",
    [
        # 2D
        (0.01, (1_000, 100), ()),
        (0.01, (1_000, 100), (None,)),
        (0.01, (1_000, 100), (slice(None),)),
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
        # 3D
        (0.005, (2_589, 101, 38), ()),
        (0.005, (2_589, 101, 38), (None,)),
        (0.005, (2_589, 101, 38), (slice(None),)),
        (0.005, (2_589, 101, 38), (slice(3, 8),)),
        (0.005, (2_589, 101, 38), (slice(3, 8), slice(45, 93))),
        (
            0.005,
            (2_589, 101, 38),
            (slice(None), slice(2, 80), [0, 22, 1, 23, 2, 24, 3]),
        ),
        # 1D
        (0.3, (10_000,), (slice(4000),)),
    ],
)
def test_blockwise_table_iter(
    a_random_sparse_nd_array: str,
    shape: Tuple[int, ...],
    coords: Tuple[Any, ...],
    a_soma_context: SOMATileDBContext,
) -> None:
    """Check blockwise iteration over non-reindexed results"""
    ndim = len(shape)
    reindex_disable_on_axis = list(range(ndim))  # disable all
    for axis, result_order in itertools.product(
        range(ndim), ["auto", "row-major", "column-major"]
    ):
        with soma.open(a_random_sparse_nd_array, mode="r", context=a_soma_context) as A:
            # get the whole enchilada in ragged form
            truth_tbl = A.read(coords=coords).tables().concat()

            block = 0
            size = max(
                5, min(A.shape[axis] // 3, 1000)
            )  # shoot for 3 blocks, in range [5,1000]
            tbls = []
            for tbl, joinids in (
                A.read(coords=coords, result_order=result_order)
                .blockwise(
                    axis=axis,
                    size=size,
                    reindex_disable_on_axis=reindex_disable_on_axis,
                )
                .tables()
            ):
                assert isinstance(tbl, pa.Table)
                assert isinstance(joinids, tuple)
                assert len(joinids) == ndim
                assert all(isinstance(joinids[d], pa.Array) for d in range(ndim))
                assert all(joinids[d].type == pa.int64() for d in range(ndim))
                assert tbl.num_columns == ndim + 1
                for d in range(ndim):
                    assert np.isin(
                        tbl.column(f"soma_dim_{d}").to_numpy(), joinids[d].to_numpy()
                    ).all()

                tbls.append(tbl)
                block += 1

            # check that stacked blocks match the full (ragged) read
            row_sort_order = [(f"soma_dim_{n}", "ascending") for n in range(ndim)]
            assert truth_tbl.sort_by(row_sort_order).equals(
                pa.concat_tables(tbls).sort_by(row_sort_order)
            )


@pytest.mark.parametrize(
    # use density to keep the sparse array memory use "reasonable"
    "density,shape",
    [
        (0.3, (10_000,)),
        (0.1, (1_000, 100)),
        (0.01, (1_000, 100, 10)),
    ],
)
@pytest.mark.parametrize("size", (999, 2**16, 2**20))
def test_blockwise_table_iter_size(
    a_random_sparse_nd_array: str, shape: Tuple[int, ...], size: int
) -> None:
    """
    Verify that blockwise iteration correctly obeys size param.
    NB: test requires soma_joinids assigned [0, n)
    """
    ndim = len(shape)
    reindex_disable_on_axis = list(range(ndim))  # reindexing off
    for axis in range(ndim):
        with soma.open(a_random_sparse_nd_array, mode="r") as A:
            assert shape == A.shape
            block = 0
            for tbl, joinids in (
                A.read()
                .blockwise(
                    axis=axis,
                    size=size,
                    reindex_disable_on_axis=reindex_disable_on_axis,
                )
                .tables()
            ):
                axis_coords = tbl.column(f"soma_dim_{axis}").to_numpy()

                assert len(joinids[axis]) <= size

                # Verify all coords are in expected range
                assert np.logical_and(
                    axis_coords >= (block * size), axis_coords < ((block + 1) * size)
                ).all()

                # Verify all block axis join ids are in same range
                assert np.logical_and(
                    joinids[axis].to_numpy() >= (block * size),
                    joinids[axis].to_numpy() < ((block + 1) * size),
                ).all()

                block += 1


@pytest.mark.parametrize(
    # use density to keep the sparse array memory use "reasonable"
    "density,shape,coords",
    [
        # 1D
        (0.01, (10_000,), ()),
        (0.01, (10_000,), (slice(200, 8000),)),
        (0.1, (10_000,), ([0, 99, 1, 100, 2, 101, 3, *list(range(150, 1000))],)),
        # 2D
        (0.0001, (1_000, 100), ()),
        (0.001, (1_000, 100), (None,)),
        (0.001, (1_000, 100), (slice(None),)),
        (0.01, (1_000, 100), (slice(10, 20),)),
        (0.01, (1_000, 100), (None, slice(10, 20))),
        (0.01, (1_000, 100), (slice(10, 20), slice(2, 33))),
        (0.01, (1_000, 100), (1,)),
        (0.01, (1_000, 100), (None, 3)),
        (0.01, (1_000, 100), (10, 3)),
        (0.01, (1_000, 100), (slice(None), 1)),
        (0.01, (1_000, 100), ([3, 99, 101, 102, 77, 0],)),
        (0.01, (1_000, 100), (slice(500), [3, 99, 33, 77, 0])),
        (0.01, (1_000, 100), (np.arange(0, 101),)),
        (0.01, (1_000, 100), (np.arange(0, 1000), [0, 10, 20])),
        (0.01, (1_000, 100), (np.arange(0, 1000), slice(2, 99))),
        # 3D
        (0.001, (1_000, 100, 10), ()),
        (0.001, (1_000, 100, 10), (slice(10, 100), slice(20, 43))),
        (0.001, (1_000, 100, 10), ([1, 2, 88, 282, 0, 382], slice(99), slice(1, 10))),
    ],
)
def test_blockwise_table_iter_reindex(
    a_random_sparse_nd_array: str,
    shape: Tuple[int, ...],
    coords: Tuple[Any, ...],
    a_soma_context: SOMATileDBContext,
) -> None:
    """Test blockwise table iteration with reindexing"""
    ndim = len(shape)
    for axis in range(ndim):
        with soma.open(a_random_sparse_nd_array, mode="r", context=a_soma_context) as A:
            # SparseCOOMatrix does not allow empty matrices.
            # See https://issues.apache.org/jira/browse/ARROW-17933
            # for more details.
            try:
                truth_coo = A.read(coords=coords).coos().concat().to_pydata_sparse()
            except pa.ArrowInvalid:
                truth_coo = None

            size = max(
                50, min(A.shape[axis] // 3, 7500)
            )  # shoot for 3 blocks, in range [50,7500]

            for tbl, joinids in (
                A.read(coords=coords).blockwise(axis=axis, size=size).tables()
            ):
                # Verify that reindex columns are [0, n)
                for d in range(ndim):
                    assert np.isin(
                        tbl.column(f"soma_dim_{d}").to_numpy(),
                        np.arange(0, len(joinids[d])),
                    ).all()

                # Check for value match with slice of whole truth. If our slice is
                # empty, skip due to aforementioned Arrow bug.
                block_coo = None
                if len(tbl) > 0:
                    d = tbl.column("soma_data").to_numpy()
                    c = np.array(
                        [tbl.column(f"soma_dim_{n}").to_numpy() for n in range(ndim)]
                    ).T
                    shape = [len(a) for a in joinids]
                    block_coo = pa.SparseCOOTensor.from_numpy(
                        d, c, shape=shape
                    ).to_pydata_sparse()

                if truth_coo is not None and block_coo is not None:
                    truth_slice = truth_coo
                    # Sparse does not allow indexing multiple dimensions with arrays, so
                    # do a progress slicing to achieve same effect.
                    for d in range(ndim):
                        truth_slice = operator.getitem(
                            truth_slice,
                            tuple([slice(None)] * d + [joinids[d].to_numpy()]),
                        )
                    assert (truth_slice == block_coo).all()


@pytest.mark.parametrize("density,shape", [(0.1, (100, 100))])
def test_blockwise_table_iter_error_checks(
    a_random_sparse_nd_array: str, shape: Tuple[int, ...]
) -> None:
    with soma.open(a_random_sparse_nd_array, mode="r") as A:
        with pytest.raises(NotImplementedError):
            next(A.read().blockwise(axis=0).tables().concat())


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
        (0.0005, (10_101, 39), ()),
        (0.0005, (10_101, 39), (None,)),
        (0.0005, (10_101, 39), (slice(None),)),
        (0.0005, (10_101, 389), (slice(33, 1048),)),
        (0.0005, (10_101, 389), (slice(33, 1048), slice(45, 333))),
    ],
)
@pytest.mark.parametrize("size", [777, 1001, 2**16])
def test_blockwise_scipy_iter(
    a_random_sparse_nd_array: str,
    coords: Tuple[Any, ...],
    size: int,
    a_soma_context: SOMATileDBContext,
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
        with soma.open(a_random_sparse_nd_array, mode="r", context=a_soma_context) as A:
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

            # Reindexing is on by default. Disable if we don't want it for minor axis.
            reindex_disable_on_axis = [minor_axis] if not reindex_sparse_axis else None
            results = []
            for sp, joinids in (
                A.read(coords)
                .blockwise(
                    axis=axis,
                    size=size,
                    reindex_disable_on_axis=reindex_disable_on_axis,
                )
                .scipy(compress=compress)
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
                assert len(joinids[axis]) == min(size, sp.shape[axis])
                assert (
                    not reindex_sparse_axis
                    or len(joinids[minor_axis]) == sp.shape[minor_axis]
                )

                # internal layout (dups, ordering, etc)
                if compress:
                    assert sp.has_canonical_format
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
def test_blockwise_scipy_iter_error_checks(
    a_random_sparse_nd_array: str, shape: Tuple[int, ...]
) -> None:
    with soma.open(a_random_sparse_nd_array, mode="r") as A:
        with pytest.raises(ValueError):
            next(A.read().blockwise(axis=2).scipy())

        with pytest.raises(soma.SOMAError):
            next(A.read().blockwise(axis=0, reindex_disable_on_axis=[0]).scipy())

        with pytest.raises(soma.SOMAError):
            next(A.read().blockwise(axis=1, reindex_disable_on_axis=[1]).scipy())


@pytest.mark.parametrize("density,shape", [(0.1, (4, 8, 16))])
def test_blockwise_scipy_iter_not_2D(
    a_random_sparse_nd_array: str, shape: Tuple[int, ...]
) -> None:
    with soma.open(a_random_sparse_nd_array, mode="r") as A:
        with pytest.raises(soma.SOMAError):
            next(A.read().blockwise(axis=0).scipy())


@pytest.mark.parametrize("density,shape", [(0.01, (10_000, 1230))])
def test_blockwise_scipy_iter_eager(
    a_random_sparse_nd_array: str,
    shape: Tuple[int, ...],
    a_soma_context: SOMATileDBContext,
) -> None:
    """Should get same results with any eager setting"""
    coords = (slice(3, 9993), slice(21, 1111))
    with soma.open(a_random_sparse_nd_array, mode="r", context=a_soma_context) as A:
        sp1 = sparse.vstack(
            [
                sp
                for sp, _ in A.read(coords)
                .blockwise(axis=0, size=1000, eager=True)
                .scipy()
            ]
        )
        sp2 = sparse.vstack(
            [
                sp
                for sp, _ in A.read(coords)
                .blockwise(axis=0, size=1000, eager=False)
                .scipy()
            ]
        )

        assert (sp1 != sp2).nnz == 0


@pytest.mark.parametrize("density,shape", [(0.001, (9799, 1530))])
def test_blockwise_scipy_iter_result_order(a_random_sparse_nd_array: str) -> None:
    """
    Confirm behavior with different result_order.
    """
    coords = (slice(7, 8693), slice(21, 999))

    with soma.open(a_random_sparse_nd_array, mode="r") as A:
        for result_order in ["auto", "row-major", "column-major"]:
            for axis in (0, 1):
                for compress in (True, False):
                    sp, _ = next(
                        A.read(coords, result_order=result_order)
                        .blockwise(axis=axis)
                        .scipy(compress=compress)
                    )

                    if compress:
                        # CS{C,R} is always sorted, always canonical
                        assert not isinstance(sp, sparse.coo_matrix)
                        assert sp.has_sorted_indices
                        assert sp.has_canonical_format
                        sp.check_format()  # raises if malformed

                    else:
                        # always canonical if row-major, regardless of format
                        if result_order == "row-major":
                            assert sp.has_canonical_format
                        assert isinstance(sp, sparse.coo_matrix)


@pytest.mark.parametrize("density,shape", [(0.001, (9799, 1530))])
@pytest.mark.parametrize(
    "coords,expected_indices",
    [
        ((), (np.arange(0, 9799), np.arange(0, 1530))),
        ((slice(43, 9000),), (np.arange(43, 9001), np.arange(0, 1530))),
        (
            (slice(430, 7100), slice(100, 200)),
            (np.arange(430, 7101), np.arange(100, 201)),
        ),
        (
            ([0, 1, 99, 3, 100, 4, 101, 5],),
            (np.array([0, 1, 99, 3, 100, 4, 101, 5]), np.arange(0, 1530)),
        ),
    ],
)
def test_blockwise_indices(
    a_random_sparse_nd_array: str,
    coords: Tuple[Any, ...],
    expected_indices: Tuple[Any, ...],
) -> None:
    """Verify indices look reasonable"""
    size = 1111

    # blockwise table
    with soma.open(a_random_sparse_nd_array, mode="r") as A:
        for axis, reindex_disable_on_axis in itertools.product(
            (0, 1), (None, [0], [1], [0, 1])
        ):
            minor_axis = 1 - axis
            block = 0
            for _, indices in (
                A.read(coords)
                .blockwise(
                    axis=axis,
                    size=size,
                    reindex_disable_on_axis=reindex_disable_on_axis,
                    eager=False,
                )
                .tables()
            ):
                assert len(indices) == 2
                beg = block * size
                assert np.array_equal(
                    expected_indices[axis][beg : beg + size], indices[axis].to_numpy()
                )
                assert np.array_equal(
                    expected_indices[minor_axis], indices[minor_axis].to_numpy()
                )
                block += 1

        # blockwise scipy
        for axis in (0, 1):
            minor_axis = 1 - axis
            for reindex_disable_on_axis in (None, [minor_axis]):
                block = 0
                for _, indices in (
                    A.read(coords)
                    .blockwise(
                        axis=axis,
                        size=size,
                        reindex_disable_on_axis=reindex_disable_on_axis,
                        eager=True,
                    )
                    .scipy()
                ):
                    assert len(indices) == 2
                    beg = block * size
                    assert np.array_equal(
                        expected_indices[axis][beg : beg + size], indices[axis]
                    )
                    assert np.array_equal(
                        expected_indices[minor_axis], indices[minor_axis]
                    )
                    block += 1


@pytest.mark.parametrize("density,shape", [(0.1, (100, 100))])
@pytest.mark.parametrize("coords", [(slice(0, 10),), (slice(1, 10),)])
def test_blockwise_scipy_reindex_disable_major_dim(
    a_random_sparse_nd_array: str, coords: Tuple[Any, ...]
) -> None:
    """
    Disable reindexing on major axis. Expected behavior:
    * fails if compress==True
    * succeeds if compress==False
    """

    with soma.open(a_random_sparse_nd_array) as A:
        for axis in (0, 1):
            # Should fail if compress==True (CSR/CSC)
            with pytest.raises(soma.SOMAError):
                next(
                    A.read(coords)
                    .blockwise(axis=axis, reindex_disable_on_axis=axis)
                    .scipy(compress=True)
                )

            # should succeed if compress==False (COO)
            sp, _ = next(
                A.read(coords)
                .blockwise(axis=axis, reindex_disable_on_axis=axis)
                .scipy(compress=False)
            )
            assert isinstance(sp, sparse.coo_matrix)


@pytest.mark.parametrize("density,shape", [(0.1, (100, 100))])
def test_blockwise_iterator_uses_thread_pool_from_context(
    a_random_sparse_nd_array: str, shape: Tuple[int, ...]
) -> None:
    pool = mock.Mock(wraps=futures.ThreadPoolExecutor(max_workers=2))
    pool.submit.assert_not_called()

    context = SOMATileDBContext(threadpool=pool)
    with soma.open(a_random_sparse_nd_array, mode="r", context=context) as A:
        axis = 0
        size = 50
        tbls = (
            A.read()
            .blockwise(
                axis=axis,
                size=size,
            )
            .tables()
        )

        # The iteration needs to happen to ensure the threadpool is used
        for tbl in tbls:
            assert tbl is not None

        pool.submit.assert_called()

    pool.reset_mock()
    pool.submit.assert_not_called()

    with soma.open(a_random_sparse_nd_array, mode="r", context=context) as A:
        axis = 0
        size = 50
        arrs = (
            A.read()
            .blockwise(
                axis=axis,
                size=size,
            )
            .scipy()
        )

        # The iteration needs to happen to ensure the threadpool is used
        for arr in arrs:
            assert arr is not None

        pool.submit.assert_called()

    pool.shutdown()


def test_global_writes(tmp_path):
    write_options = soma.TileDBWriteOptions(**{"sort_coords": False})

    soma.SparseNDArray.create(tmp_path.as_posix(), type=pa.uint8(), shape=(3,))

    with pytest.raises(
        soma.SOMAError,
        match=r"Write failed; Coordinates (.*) succeed (.*) in the global order",
    ):
        with soma.SparseNDArray.open(tmp_path.as_posix(), "w") as A:
            A.write(
                pa.Table.from_pydict(
                    {
                        "soma_dim_0": pa.array([2, 1, 0], type=pa.int64()),
                        "soma_data": pa.array([1, 2, 3], type=pa.uint8()),
                    }
                ),
                platform_config=write_options,
            )

    data = pa.Table.from_pydict(
        {
            "soma_dim_0": pa.array([0, 1, 2], type=pa.int64()),
            "soma_data": pa.array([1, 2, 3], type=pa.uint8()),
        }
    )

    with soma.SparseNDArray.open(tmp_path.as_posix(), "w") as A:
        A.write(
            data,
            platform_config=write_options,
        )

    with soma.SparseNDArray.open(tmp_path.as_posix()) as A:
        assert A.read().tables().concat() == data

    with pytest.raises(ValueError):
        # Takes TileDBWriteOptions as of TileDB-SOMA 1.13
        with soma.SparseNDArray.open(tmp_path.as_posix(), "w") as A:
            A.write(
                data,
                platform_config=soma.TileDBCreateOptions(),
            )


def test_pass_configs(tmp_path):
    uri = tmp_path.as_posix()

    with soma.SparseNDArray.create(
        tmp_path.as_posix(), type=pa.uint8(), shape=(3,)
    ) as a:
        data = pa.Table.from_pydict(
            {
                "soma_dim_0": pa.array([0, 1, 2], type=pa.int64()),
                "soma_data": pa.array([1, 2, 3], type=pa.uint8()),
            }
        )
        a.write(data)

    # Pass a custom config to open
    with soma.SparseNDArray.open(
        uri,
        "r",
        context=soma.SOMATileDBContext(
            {"sm.mem.total_budget": "0", "sm.io_concurrency_level": "0"}
        ),
    ) as sdf:

        # This errors out as 0 is not a valid value to set the total memory
        # budget or number of threads
        with pytest.raises(soma.SOMAError):
            next(sdf.read().tables())

        # This still errors out because read still sees that the number of
        # threads is 0 and therefore invalid
        with pytest.raises(soma.SOMAError):
            next(sdf.read(platform_config={"sm.mem.total_budget": "10000"}).tables())

        # With correct values, this reads without issue
        next(
            sdf.read(
                platform_config={
                    "sm.mem.total_budget": "10000",
                    "sm.io_concurrency_level": "1",
                }
            ).tables()
        )


def test_iter(tmp_path: pathlib.Path):
    arrow_tensor = create_random_tensor("table", (1,), np.float32(), density=1)

    with soma.SparseNDArray.create(
        tmp_path.as_uri(), type=pa.float64(), shape=(1,)
    ) as write_arr:
        write_arr.write(arrow_tensor)

    # Verify that the SOMAArray stays open as long as the ManagedQuery
    # (i.e., `next`) is still active
    a = soma.open(tmp_path.as_uri(), mode="r").read().tables()
    assert next(a)
    with pytest.raises(StopIteration):
        next(a)

    # Open two instances of the same array. Iterating through one should not
    # affect the other
    a = soma.open(tmp_path.as_uri(), mode="r").read().tables()
    b = soma.open(tmp_path.as_uri(), mode="r").read().tables()
    assert next(a)
    assert next(b)
    with pytest.raises(StopIteration):
        next(a)
    with pytest.raises(StopIteration):
        next(b)


def test_context_cleanup(tmp_path: pathlib.Path) -> None:
    arrow_tensor = create_random_tensor("table", (1,), np.float32(), density=1)
    with soma.SparseNDArray.create(
        tmp_path.as_uri(), type=pa.float64(), shape=(1,)
    ) as write_arr:
        write_arr.write(arrow_tensor)

    def test(path, tiledb_config):
        context = soma.SOMATileDBContext().replace(tiledb_config=tiledb_config)
        X = soma.SparseNDArray.open(path, context=context, mode="r")
        mq = X.read().tables()
        return mq

    for _ in range(100):
        # Run test multiple times. While the C++ this tests (dtor order)
        # is deterministic, it is triggered by the Python GC, which makes
        # no specific guarantees about when it will sweep any given object.
        _ = test(
            tmp_path.as_uri(),
            {
                "vfs.s3.no_sign_request": "true",
                "vfs.s3.region": "us-west-2",
            },
        )
        gc.collect()


def test_sparse_nd_array_null(tmp_path):
    uri = tmp_path.as_posix()

    pydict = {
        "soma_dim_0": pa.array([None, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
        "soma_data": pa.array(
            [None, 0, None, 1, 2, None, None, 3, 4, 5], type=pa.float64()
        ),
    }
    table = pa.Table.from_pydict(pydict)

    soma.SparseNDArray.create(uri, type=pa.int64(), shape=(10,))

    # As of version 1.15.6 we were throwing in this case. However, we found
    # a compatibility issue with pyarrow versions below 17. Thus this is
    # now non-fatal.
    # with soma.SparseNDArray.open(uri, "w") as A:
    #    with raises_no_typeguard(soma.SOMAError):
    #        # soma_joinid cannot be nullable
    #        A.write(table[:5])
    #        A.write(table[5:])

    pydict["soma_dim_0"] = pa.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    table = pa.Table.from_pydict(pydict)

    with soma.SparseNDArray.open(uri, "w") as A:
        A.write(table[:5])
        A.write(table[5:])

    with soma.SparseNDArray.open(uri) as A:
        pdf = A.read().tables().concat()

        # soma_data is a non-nullable attribute. In ManagedQuery.set_array_data,
        # any null values present in non-nullable attributes get casted to
        # fill values. In the case for float64, the fill value is 0
        np.testing.assert_array_equal(pdf["soma_data"], table["soma_data"].fill_null(0))


def test_reopen_metadata_sc61118(tmp_path):
    uri = tmp_path.as_posix()
    with soma.SparseNDArray.create(uri, type=pa.int64(), shape=(10,)) as A1:
        A1.metadata["foo"] = "bar"
        with A1.reopen(mode="r") as A2:
            assert dict(A1.metadata) == dict(A2.metadata)


def test_reopen_shape_sc61123(tmp_path):
    uri = tmp_path.as_posix()
    with soma.SparseNDArray.create(uri, type=pa.int64(), shape=(10,)) as A:
        assert A.shape == (10,)
        assert isinstance(A, soma.SparseNDArray)
        A = A.reopen(mode="r")
        assert A.shape == (10,)
        assert isinstance(A, soma.SparseNDArray)
