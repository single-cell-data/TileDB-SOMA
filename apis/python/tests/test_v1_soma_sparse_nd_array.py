# import pytest

import numpy as np
import pyarrow as pa

import tiledbsc.v1


def test_soma_sparse_nd_array_ok_no_storage():
    arr = tiledbsc.v1.SOMASparseNdArray(uri="/foo/bar")
    assert arr.get_uri() == "/foo/bar"
    assert arr.get_name() == "bar"


def test_soma_sparse_nd_array(tmp_path):
    arr = tiledbsc.v1.SOMASparseNdArray(uri=tmp_path.as_posix())

    nr = 10
    nc = 20
    arr.create(pa.float64(), [nr, nc])

    tensor = pa.SparseCOOTensor.from_numpy(
        data=np.asarray([7.0, 8.0, 9.0]),
        coords=[[1, 2], [3, 4], [5, 6]],
        shape=(nr, nc),
    )
    arr.write(tensor)

    # TODO: check more things
