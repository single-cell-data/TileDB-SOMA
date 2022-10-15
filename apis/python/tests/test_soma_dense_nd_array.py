from typing import Tuple

import numpy as np
import pyarrow as pa
import pytest

import tiledbsoma as soma

from . import NDARRAY_ARROW_TYPES_NOT_SUPPORTED, NDARRAY_ARROW_TYPES_SUPPORTED


def test_soma_dense_nd_array_ok_no_storage():
    arr = soma.SOMADenseNdArray(uri="/foo/bar/")
    assert arr.uri == "/foo/bar/"
    assert not arr.exists()
    assert arr.soma_type == "SOMADenseNdArray"


@pytest.mark.parametrize(
    "shape", [(10,), (1, 100), (10, 1, 100), (2, 4, 6, 8), [1], (1, 2, 3, 4, 5)]
)
@pytest.mark.parametrize("element_type", NDARRAY_ARROW_TYPES_SUPPORTED)
def test_soma_dense_nd_array_create_ok(
    tmp_path, shape: Tuple[int, ...], element_type: pa.DataType
):
    """
    Test all cases we expect "create" to succeed.
    """
    assert pa.types.is_primitive(element_type)  # sanity check incoming params

    a = soma.SOMADenseNdArray(uri=tmp_path.as_posix())
    a.create(element_type, shape)
    assert a.soma_type == "SOMADenseNdArray"
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
def test_soma_dense_nd_array_create_fail(
    tmp_path, shape: Tuple[int, ...], element_type: pa.DataType
):
    a = soma.SOMADenseNdArray(uri=tmp_path.as_posix())
    with pytest.raises(TypeError):
        a.create(element_type, shape)
    assert not a.exists()


def test_soma_dense_nd_array_delete(tmp_path):
    a = soma.SOMADenseNdArray(uri=tmp_path.as_posix())
    a.create(pa.int8(), (100, 100))
    assert a.exists()

    a.delete()
    assert not a.exists()

    # should be silent about non-existent object
    assert a.delete() is None
    assert soma.SOMADenseNdArray(uri="no such array").delete() is None


@pytest.mark.parametrize("shape", [(10,), (10, 20), (10, 20, 2), (2, 4, 6, 8)])
def test_soma_dense_nd_array_read_write_tensor(tmp_path, shape: Tuple[int, ...]):
    a = soma.SOMADenseNdArray(tmp_path.as_posix())
    a.create(pa.float64(), shape)
    ndim = len(shape)

    # random sample - written to entire array
    data = np.random.default_rng().standard_normal(np.prod(shape)).reshape(shape)
    coords = tuple(slice(0, dim_len) for dim_len in shape)
    a.write_tensor(coords, pa.Tensor.from_numpy(data))
    del a

    # check multiple read paths
    b = soma.SOMADenseNdArray(tmp_path.as_posix())

    t = b.read_tensor((slice(None),) * ndim, result_order="row-major")
    assert t.equals(pa.Tensor.from_numpy(data))

    t = b.read_tensor((slice(None),) * ndim, result_order="column-major")
    assert t.equals(pa.Tensor.from_numpy(data.transpose()))

    # write a single-value sub-array and recheck
    b.write_tensor(
        (0,) * len(shape),
        pa.Tensor.from_numpy(np.zeros((1,) * len(shape), dtype=np.float64)),
    )
    data[(0,) * len(shape)] = 0.0
    t = b.read_tensor((slice(None),) * ndim)
    assert t.equals(pa.Tensor.from_numpy(data))


@pytest.mark.parametrize("shape", [(), (0,), (10, 0), (0, 10), (1, 2, 0)])
def test_zero_length_fail(tmp_path, shape):
    """Zero length dimensions are expected to fail"""
    a = soma.SOMADenseNdArray(tmp_path.as_posix())
    with pytest.raises(ValueError):
        a.create(type=pa.float32(), shape=shape)


def test_soma_dense_nd_array_reshape(tmp_path):
    """
    Reshape currently unimplemented.
    """
    a = soma.SOMADenseNdArray(tmp_path.as_posix())
    a.create(type=pa.int32(), shape=(10, 10, 10))
    with pytest.raises(NotImplementedError):
        assert a.reshape((100, 10, 1))
