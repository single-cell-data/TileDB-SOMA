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
    def np_dense_nd_array(self):
        arr = soma.DenseNDArray.create(self.uri, type=pa.int32(), shape=(2, 3, 4))
        ndarray = np.random.default_rng().integers(0, 10, 24).reshape(2, 3, 4)
        data = pa.Tensor.from_numpy(ndarray)
        arr.write((slice(None),), data)
        return ndarray

    def base_R_script(self):
        return f"""
        library("tiledbsoma")
        soma_ndarray <- SOMADenseNDArrayOpen("{self.uri}")
        table <- soma_ndarray$read_arrow_table()
        df <- as.data.frame(table)
        """

    def test_ndarray_shape_matches(self, np_dense_nd_array):
        """
        The source ndarray is a 2x3x4 tensor, so the resulting soma_ndarray should have 3 dimensions.
        """
        self.r_assert("stopifnot(length(soma_ndarray$dimensions()) == 3)")

    def test_ndarray_type_matches(self, np_dense_nd_array):
        """
        The DenseNDArray should have a type of int32.
        """
        self.r_assert('stopifnot(table$soma_data$type$ToString() == "int32")')

    def test_ndarray_content_matches(self, np_dense_nd_array):
        """
        The DenseNDArray content should match. Sparse reads materialize coordinates, e.g.
        soma_dim_0, soma_dim_1, soma_dim_2, and soma_data. Dense reads do not: we get soma_data.
        A quick way to check for matches is to flatten out both sides and compare indexwise.
        """
        flattened_python_read = np_dense_nd_array.flatten().tolist()
        multi_assert = ""
        for i in range(2 * 3 * 4):
            val = flattened_python_read[i]
            multi_assert += f"stopifnot(df$soma_data[[{i}+1]] == {val})"
            multi_assert += "\n"
        self.r_assert(multi_assert)
