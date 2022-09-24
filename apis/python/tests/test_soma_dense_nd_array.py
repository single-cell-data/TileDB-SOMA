from typing import Tuple

import numpy as np
import pyarrow as pa
import pytest

import tiledbsoma as soma

from . import NDARRAY_ARROW_TYPES_NOT_SUPPORTED, NDARRAY_ARROW_TYPES_SUPPORTED


"""
TODO:
- [X] create
- [-] delete
- [X] exists
- [X] get type
- [X] get shape
- [X] get ndims
- [X] get schema
- [X] get is_sparse
- [ ] get metadata
- [ ] read
- [ ] write
- [ ] reshape
"""


def test_soma_dense_nd_array_ok_no_storage():
    arr = soma.SOMADenseNdArray(uri="/foo/bar")
    assert arr.uri == "/foo/bar"
    assert arr.name == "bar"
    assert not arr.exists()
    assert arr.type == "SOMADenseNdArray"


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
    assert a.type == "SOMADenseNdArray"
    assert a.uri == tmp_path.as_posix()
    assert a.ndims == len(shape)
    assert a.shape == tuple(shape)
    assert a.is_sparse is False
    assert a.exists()

    assert a.schema is not None
    expected_field_names = ["data"] + [f"__dim_{d}" for d in range(len(shape))]
    assert set(a.schema.names) == set(expected_field_names)
    for d in range(len(shape)):
        assert a.schema.field(f"__dim_{d}").type == pa.uint64()
    assert a.schema.field("data").type == element_type


@pytest.mark.parametrize("shape", [(10,)])
@pytest.mark.parametrize("element_type", NDARRAY_ARROW_TYPES_NOT_SUPPORTED)
def test_soma_dense_nd_array_create_fail(
    tmp_path, shape: Tuple[int, ...], element_type: pa.DataType
):
    a = soma.SOMADenseNdArray(uri=tmp_path.as_posix())
    with pytest.raises(TypeError):
        a.create(element_type, shape)
    assert not a.exists()


# TODO: implement test when operation is implemented
@pytest.mark.xfail(reason="Unimplemented operation")
def test_soma_dense_nd_array_delete(tmp_path):
    a = soma.SOMADenseNdArray(uri=tmp_path.as_posix())
    a.create(pa.int8(), (100, 100))
    assert a.exists()

    a.delete()
    assert not a.exists()




""" XXX obsolete - remove """
def test_soma_dense_nd_array(tmp_path):
    nr = 10
    nc = 20
    a = soma.SOMADenseNdArray(tmp_path.as_posix())

    a.create(pa.float64(), [nr, nc])

    a.write((slice(0, nr), slice(0, nc)), pa.Tensor.from_numpy(np.eye(nr, nc)))
    # a.write((slice(8, 12), slice(10, 16)), pa.Tensor.from_numpy(np.ones((4, 6))))

    # TODO: check more things
