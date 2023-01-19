from typing import Tuple

import numpy as np
import pyarrow as pa
import pytest
import tiledb

import tiledbsoma as soma
from tiledbsoma.options import SOMATileDBContext

from . import NDARRAY_ARROW_TYPES_NOT_SUPPORTED, NDARRAY_ARROW_TYPES_SUPPORTED


def test_dense_nd_array_ok_no_storage():
    arr = soma.DenseNDArray(uri="/foo/bar/")
    assert arr.uri == "/foo/bar/"
    assert not arr.exists()
    assert arr.soma_type == "SOMADenseNDArray"


@pytest.mark.parametrize(
    "shape", [(10,), (1, 100), (10, 1, 100), (2, 4, 6, 8), [1], (1, 2, 3, 4, 5)]
)
@pytest.mark.parametrize("element_type", NDARRAY_ARROW_TYPES_SUPPORTED)
def test_dense_nd_array_create_ok(
    tmp_path, shape: Tuple[int, ...], element_type: pa.DataType
):
    """
    Test all cases we expect "create" to succeed.
    """
    assert pa.types.is_primitive(element_type)  # sanity check incoming params

    a = soma.DenseNDArray(uri=tmp_path.as_posix())
    a.create(element_type, shape)
    assert a.soma_type == "SOMADenseNDArray"
    assert a.uri == tmp_path.as_posix()
    assert a.ndim == len(shape)
    assert a.shape == tuple(shape)
    assert a.is_sparse is False
    assert a.exists()

    assert a.schema is not None
    expected_field_names = ["soma_data"] + [f"soma_dim_{d}" for d in range(len(shape))]
    assert set(a.schema.names) == set(expected_field_names)
    for d in range(len(shape)):
        assert a.schema.field(f"soma_dim_{d}").type == pa.int64()
    assert a.schema.field("soma_data").type == element_type


@pytest.mark.parametrize("shape", [(10,)])
@pytest.mark.parametrize("element_type", NDARRAY_ARROW_TYPES_NOT_SUPPORTED)
def test_dense_nd_array_create_fail(
    tmp_path, shape: Tuple[int, ...], element_type: pa.DataType
):
    a = soma.DenseNDArray(uri=tmp_path.as_posix())
    with pytest.raises(TypeError):
        a.create(element_type, shape)
    assert not a.exists()


def test_dense_nd_array_delete(tmp_path):
    a = soma.DenseNDArray(uri=tmp_path.as_posix())
    a.create(pa.int8(), (100, 100))
    assert a.exists()

    a.delete()
    assert not a.exists()

    # should be silent about non-existent object
    assert a.delete() is None
    assert soma.DenseNDArray(uri="no such array").delete() is None


@pytest.mark.parametrize("shape", [(10,), (10, 20), (10, 20, 2), (2, 4, 6, 8)])
def test_dense_nd_array_read_write_tensor(tmp_path, shape: Tuple[int, ...]):
    a = soma.DenseNDArray(tmp_path.as_posix())
    a.create(pa.float64(), shape)
    ndim = len(shape)

    # random sample - written to entire array
    data = np.random.default_rng().standard_normal(np.prod(shape)).reshape(shape)
    coords = tuple(slice(0, dim_len) for dim_len in shape)
    a.write(coords, pa.Tensor.from_numpy(data))
    del a

    # check multiple read paths
    b = soma.DenseNDArray(tmp_path.as_posix())

    t = b.read((slice(None),) * ndim, result_order="row-major")
    assert t.equals(pa.Tensor.from_numpy(data))

    t = b.read((slice(None),) * ndim, result_order="column-major")
    assert t.equals(pa.Tensor.from_numpy(data.transpose()))

    # write a single-value sub-array and recheck
    b.write(
        (0,) * len(shape),
        pa.Tensor.from_numpy(np.zeros((1,) * len(shape), dtype=np.float64)),
    )
    data[(0,) * len(shape)] = 0.0
    t = b.read((slice(None),) * ndim)
    assert t.equals(pa.Tensor.from_numpy(data))


@pytest.mark.parametrize("shape", [(), (0,), (10, 0), (0, 10), (1, 2, 0)])
def test_zero_length_fail(tmp_path, shape):
    """Zero length dimensions are expected to fail"""
    a = soma.DenseNDArray(tmp_path.as_posix())
    with pytest.raises(ValueError):
        a.create(type=pa.float32(), shape=shape)


def test_dense_nd_array_reshape(tmp_path):
    """
    Reshape currently unimplemented.
    """
    a = soma.DenseNDArray(tmp_path.as_posix())
    a.create(type=pa.int32(), shape=(10, 10, 10))
    with pytest.raises(NotImplementedError):
        assert a.reshape((100, 10, 1))


@pytest.mark.parametrize(
    "io",
    [
        {
            "coords": (2, 3),
            "output": np.array([[203]]),
        },
        {
            "coords": (slice(None), 3),
            "output": np.array([[3], [103], [203], [303]]),
        },
        {
            "coords": (2, slice(None)),
            "output": np.array([[200, 201, 202, 203, 204, 205]]),
        },
        {
            "coords": (slice(None, 2), slice(5, None)),
            "output": np.array([[5], [105], [205]]),
        },
        {
            "coords": (slice(0, 2), slice(5, 5)),
            "output": np.array([[5], [105], [205]]),
        },
        {
            "coords": (slice(None), slice(None)),
            "cfg": {
                "soma.init_buffer_bytes": 100
            },  # Known small enough to force multiple reads
            "output": np.array(
                [
                    [0, 1, 2, 3, 4, 5],
                    [100, 101, 102, 103, 104, 105],
                    [200, 201, 202, 203, 204, 205],
                    [300, 301, 302, 303, 304, 305],
                ]
            ),
        },
    ],
)
def test_dense_nd_array_slicing(tmp_path, io):
    """
    We already have tests that check n-d for various values of n. This one (which happens to use 2-d
    data, though not in an essential way) checks subarray slicing. In particular, it validates
    SOMA's doubly-inclusive slice indexing semantics against Python's singly-inclusive slicing
    semantics, ensuring that none of the latter has crept into the former.
    """
    cfg = {}
    if "cfg" in io:
        cfg = io["cfg"]
    context = SOMATileDBContext(tiledb_ctx=tiledb.Ctx(cfg))

    a = soma.DenseNDArray(tmp_path.as_posix(), context=context)
    nr = 4
    nc = 6

    a.create(pa.int64(), [nr, nc])
    npa = np.zeros((nr, nc))
    for i in range(nr):
        for j in range(nc):
            npa[i, j] = 100 * i + j
    a.write(coords=(slice(0, nr), slice(0, nc)), values=pa.Tensor.from_numpy(npa))

    if "throws" in io:
        with pytest.raises(io["throws"]):
            a.read(io["coords"]).to_numpy()
    else:
        output = a.read(io["coords"]).to_numpy()
        assert np.all(output == io["output"])


@pytest.mark.parametrize(
    "io",
    [
        {
            "shape": (10,),
            "coords": (-1,),
            "throws": (RuntimeError, tiledb.cc.TileDBError),
        },
        {
            "shape": (10,),
            "coords": (12,),
            "throws": (RuntimeError, tiledb.cc.TileDBError),
        },
        {
            "shape": (10,),
            "coords": (
                2,
                3,
            ),
            "throws": ValueError,
        },
        {
            "shape": (10, 20),
            "coords": (
                2,
                3,
                4,
            ),
            "throws": ValueError,
        },
        {
            "shape": (10, 20),
            "coords": (slice(-2, -1),),
            "throws": ValueError,
        },
        {
            "shape": (10, 20),
            "coords": (slice(2, 3, -1),),
            "throws": ValueError,
        },
        {
            "shape": (10, 20),
            "coords": (slice(3, 2, 1),),
            "throws": ValueError,
        },
        {
            "shape": (10, 20),
            "coords": (slice(4, 8, 2),),
            "throws": ValueError,
        },
        pytest.param(  # pending https://github.com/single-cell-data/TileDB-SOMA/issues/584
            {
                "shape": (10, 20),
                "coords": (
                    pa.array(
                        [1, 2, 3],
                    )
                ),
                "throws": ValueError,
            },
            marks=pytest.mark.xfail,
        ),
    ],
)
def test_dense_nd_array_indexing_errors(tmp_path, io):
    shape = io["shape"]
    read_coords = io["coords"]

    a = soma.DenseNDArray(tmp_path.as_posix())
    a.create(pa.int64(), shape)

    npa = np.random.default_rng().standard_normal(np.prod(shape)).reshape(shape)

    write_coords = tuple(slice(0, dim_len) for dim_len in shape)
    a.write(coords=write_coords, values=pa.Tensor.from_numpy(npa))

    with pytest.raises(io["throws"]):
        a.read(coords=read_coords).to_numpy()
