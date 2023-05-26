import pyarrow as pa
import pytest
from scipy import sparse

import tiledbsoma as soma

from .common import TestWritePythonReadR


class TestSparseNDArrayWritePythonReadR(TestWritePythonReadR):
    """
    Tests that a SOMASparseNDArray can be written by Python and read by R.
    """

    @pytest.fixture(scope="class")
    def sparse_nd_array(self):
        arr = soma.SparseNDArray.create(self.uri, type=pa.int32(), shape=(5, 5))
        X = sparse.diags([1, 1, 1, 1], offsets=1).tocoo()
        tensor = pa.SparseCOOTensor.from_scipy(X)
        arr.write(tensor)
        return X

    def base_R_script(self):
        return f"""
        library("tiledbsoma")
        soma_ndarray <- SOMASparseNDArrayOpen("{self.uri}")
        table <- soma_ndarray$read()$tables()$concat()
        M <-  soma_ndarray$read()$sparse_matrix(zero_based=T)$concat()$get_one_based_matrix()
        df <- as.data.frame(table)
        """

    def test_ndarray_shape_matches(self, sparse_nd_array):
        """
        The source ndarray is a 5x5 sparse matrix, so the resulting soma_ndarray should have 2 dimensions.
        """
        self.r_assert("stopifnot(length(soma_ndarray$dimensions()) == 2)")
        self.r_assert("stopifnot(M@Dim == c(5, 5))")

    def test_ndarray_type_matches(self, sparse_nd_array):
        """
        The SparseNDArray should have a type of int32.
        """
        self.r_assert('stopifnot(table$soma_data$type$ToString() == "int32")')

    def test_ndarray_content_matches(self, sparse_nd_array):
        """
        The SparseNDArray content should match.
        """
        # Original matrix:
        # [0., 1., 0., 0., 0.],
        # [0., 0., 1., 0., 0.],
        # [0., 0., 0., 1., 0.],
        # [0., 0., 0., 0., 1.],
        # [0., 0., 0., 0., 0.]

        self.r_assert(
            "stopifnot(all.equal(as.matrix(M), matrix(c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0), nrow=5, ncol=5, byrow=TRUE)))"
        )
