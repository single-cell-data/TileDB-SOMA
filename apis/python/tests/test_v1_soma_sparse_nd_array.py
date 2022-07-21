# import pytest

import tiledbsc.v1


def test_soma_sparse_nd_array_ok_no_storage():
    arr = tiledbsc.v1.SOMASparseNdArray(uri="/foo/bar")
    assert arr.get_uri() == "/foo/bar"
    assert arr.get_name() == "bar"


# TODO: create and then read from storage:
#    assert arr.get_shape() == (1, 2, 3)
#    assert arr.get_ndims() == 3
#    assert arr.get_is_sparse()
#
#
# def test_soma_sparse_nd_array_errors_no_storage():
#    with pytest.raises(AssertionError):
#        tiledbsc.v1.SOMASparseNdArray(uri="/foo/bar", name="bar", shape=())
#
#    with pytest.raises(AssertionError):
#        tiledbsc.v1.SOMASparseNdArray(uri="/foo/bar", name="bar", shape=(1, 0, 3))
