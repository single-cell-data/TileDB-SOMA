import numpy as np
import pytest

import tiledbsoma as soma

from .common import TestReadPythonWriteR


class TestSparseNDArrayWriteRReadPython(TestReadPythonWriteR):
    """
    Tests that a SOMASparseNDArray can be written by R and read by Python.
    """

    @pytest.fixture(scope="class")
    def R_ndarray(self):
        base_script = f"""
        library("tiledbsoma")
        library("arrow")
        library("Matrix")

        sndarray <- SOMASparseNDArrayCreate("{self.uri}", int32(), c(5, 5))

        M <- Matrix(nrow = 5, ncol = 5, data = 0, sparse = TRUE)
        M <- as(M, "TsparseMatrix")
        M[1,2] <- 1
        M[2,3] <- 2
        sndarray$write(M)
        """
        self.execute_R_script(base_script)

    def test_ndarray_shape_matches(self, R_ndarray):
        """
        The source ndarray is a 5x5 matrix, so the resulting soma_ndarray should have 2 dimensions.
        """
        with soma.open(self.uri) as sdf:
            ndarr = sdf.read().coos().concat()
            assert ndarr.shape == (5, 5)

    def test_ndarray_type_matches(self, R_ndarray):
        """
        The SparseNDArray should have a type of int32.
        """
        with soma.open(self.uri) as sdf:
            ndarr = sdf.read().coos().concat()
            assert ndarr.type == "int32"

    def test_ndarray_content_matches(self, R_ndarray):
        """
        The SparseNDArray content should match.
        """
        with soma.open(self.uri) as sdf:
            arr = sdf.read().coos().concat().to_scipy().todense()
            assert np.array_equal(arr[0], np.matrix([0, 1, 0, 0, 0]))
            assert np.array_equal(arr[1], np.matrix([0, 0, 2, 0, 0]))
            assert np.array_equal(arr[2], np.matrix([0, 0, 0, 0, 0]))
            assert np.array_equal(arr[3], np.matrix([0, 0, 0, 0, 0]))
            assert np.array_equal(arr[4], np.matrix([0, 0, 0, 0, 0]))
