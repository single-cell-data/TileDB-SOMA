import tiledbsc.v1
import pyarrow as pa
import numpy as np

def test_soma_dense_nd_array_ok_no_storage():
    arr = tiledbsc.v1.SOMASparseNdArray(uri="/foo/bar")
    assert arr.get_uri() == "/foo/bar"
    assert arr.get_name() == "bar"

def test_soma_dense_nd_array(tmp_path):
    nr = 10
    nc = 20
    a = tiledbsc.v1.SOMADenseNdArray(tmp_path.as_posix())

    a.create(pa.float64(), [nr, nc])

    a.write((slice(0, nr), slice(0, nc)), pa.Tensor.from_numpy(np.eye(nr, nc)))
    # a.write((slice(8, 12), slice(10, 16)), pa.Tensor.from_numpy(np.ones((4, 6))))

    # TODO: check more things
