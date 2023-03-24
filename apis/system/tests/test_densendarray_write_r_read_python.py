import numpy as np
import pytest

import tiledbsoma as soma

from .common import TestReadPythonWriteR


class TestDenseNDArrayWriteRReadPython(TestReadPythonWriteR):
    """
    Tests that a SOMADenseNDArray can be written by R and read by Python.
    """

    @pytest.fixture(scope="class")
    def R_ndarray(self):
        base_script = f"""
        library("tiledbsoma")
        library("arrow")

        sndarray <- SOMADenseNDArrayCreate("{self.uri}", int32(), c(3, 2))
        mat <- matrix(c(1,2,3,4,5,6), nrow=3, ncol=2, byrow=TRUE)
        sndarray$write(mat)
        """
        self.execute_R_script(base_script)

    def test_ndarray_shape_matches(self, R_ndarray):
        """
        The source ndarray is a 3x2 matrix, so the resulting soma_ndarray should have 2 dimensions.
        """
        with soma.open(self.uri) as sdf:
            ndarr = sdf.read()
            assert ndarr.shape == (3, 2)

    def test_ndarray_type_matches(self, R_ndarray):
        """
        The DenseNDArray should have a type of int32.
        """
        with soma.open(self.uri) as sdf:
            ndarr = sdf.read()
            assert ndarr.type == "int32"

    def test_ndarray_content_matches(self, R_ndarray):
        """
        The DenseNDArray content should match.
        """
        with soma.open(self.uri) as sdf:
            arr = sdf.read().to_numpy()
            assert np.array_equal(arr[0], np.asarray([1, 2]))
            assert np.array_equal(arr[1], np.asarray([3, 4]))
            assert np.array_equal(arr[2], np.asarray([5, 6]))
