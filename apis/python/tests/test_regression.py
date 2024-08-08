"""Testing module for regression tests"""

import numpy as np
import pyarrow as pa

import tiledbsoma as soma


def test_nd_dense_non_contiguous_write(tmp_path):
    """Test regression dected in GitHub Issue #2537"""
    # Create data.
    data = (
        np.arange(np.product(24), dtype=np.uint8)
        .reshape((4, 3, 2))
        .transpose((2, 0, 1))
    )
    coords = tuple(slice(0, dim_len) for dim_len in data.shape)
    tensor = pa.Tensor.from_numpy(data)

    # Create array and write data to it.
    with soma.DenseNDArray.create(
        tmp_path.as_posix(), type=pa.uint8(), shape=data.shape
    ) as array:
        array.write(coords, tensor)

    # Check the data is correct when we read it back.
    with soma.DenseNDArray.open(tmp_path.as_posix()) as array:
        result = array.read(coords)
    np.testing.assert_equal(data, result.to_numpy())

    # Check the data is correct when we read it back.
    with soma.DenseNDArray.open(tmp_path.as_posix()) as array:
        result = array.read(coords, result_order="column-major")
    np.testing.assert_equal(data.transpose(), result.to_numpy())
