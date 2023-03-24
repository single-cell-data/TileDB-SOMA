import numpy as np
import pyarrow as pa
import pytest

import tiledbsoma as soma

from .common import TestWritePythonReadR


class TestDenseNDArrayWritePythonReadR(TestWritePythonReadR):
    """
    Tests that a SOMADenseNDArray can be written by Python and read by R.
    """

    @pytest.fixture(scope="class")
    def dense_nd_array(self):
        arr = soma.DenseNDArray.create(self.uri, type=pa.int32(), shape=(2, 3, 4))
        ndarray = np.random.default_rng().integers(0, 10, 24).reshape(2, 3, 4)
        data = pa.Tensor.from_numpy(ndarray)
        arr.write((slice(None),), data)
        return ndarray

    def base_R_script(self):
        return f"""
        library("tiledbsoma")
        soma_ndarray <- SOMADenseNDArrayOpen("{self.uri}")
        table = soma_ndarray$read_arrow_table()
        df = as.data.frame(table)
        """

    def test_ndarray_shape_matches(self, dense_nd_array):
        """
        The source ndarray is a 2x3x4 tensor, so the resulting soma_ndarray should have 3 dimensions.
        """
        self.r_assert("stopifnot(length(soma_ndarray$dimensions()) == 3)")

    def test_ndarray_type_matches(self, dense_nd_array):
        """
        The DenseNDArray should have a type of int32.
        """
        self.r_assert('stopifnot(table$soma_data$type$ToString() == "int32")')

    def test_ndarray_content_matches(self, dense_nd_array):
        """
        The DenseNDArray content should match. To test this, we convert the DenseNDArray into an arrow.Table and
        to a data.frame, and we use the coordinates to match the values.
        """
        multi_assert = ""
        for x, y, z in np.ndindex((2, 3, 4)):
            val = dense_nd_array[x, y, z]
            multi_assert += f"stopifnot(df[which(df$soma_dim_0 == {x} & df$soma_dim_1 == {y} & df$soma_dim_2 == {z}),]$soma_data == {val})"
            multi_assert += "\n"
        self.r_assert(multi_assert)
