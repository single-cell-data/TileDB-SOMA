import abc
from typing import TypeVar

import numpy as np
import pyarrow as pa
import scipy.sparse as sp
import somacore

# This package's pybind11 code
import tiledbsoma.libtiledbsoma as clib

from .exception import SOMAError
from .types import NTuple


class TableReadIter(somacore.ReadIter[pa.Table]):
    """Iterator over Arrow Table elements"""

    def __init__(self, sr: clib.SOMAReader):
        self.sr = sr

    def __next__(self) -> pa.Table:
        arrow_table = self.sr.read_next()
        if arrow_table is None:
            raise StopIteration

        return arrow_table

    def concat(self) -> pa.Table:
        """Concatenate remainder of iterator, and return as a single Arrow Table"""
        return pa.concat_tables(self)


RT = TypeVar("RT")


class SparseTensorReadIterBase(somacore.ReadIter[RT], metaclass=abc.ABCMeta):
    """Private implementation class"""

    def __init__(self, sr: clib.SOMAReader, shape: NTuple):
        self.sr = sr
        self.shape = shape

    @abc.abstractmethod
    def _from_table(self, arrow_table: pa.Table) -> RT:
        raise NotImplementedError()

    def __next__(self) -> RT:
        arrow_table = self.sr.read_next()
        if arrow_table is None:
            raise StopIteration

        return self._from_table(arrow_table)

    def concat(self) -> RT:
        """Returns all the requested data in a single operation.

        If some data has already been retrieved using ``next``, this will return
        the rest of the data after that is already returned.
        """
        arrow_tables = pa.concat_tables(TableReadIter(self.sr))
        return self._from_table(arrow_tables)


class SparseCOOTensorReadIter(SparseTensorReadIterBase[pa.SparseCOOTensor]):
    """Iterator over Arrow SparseCOOTensor elements"""

    def _from_table(self, arrow_table: pa.Table) -> pa.SparseCOOTensor:
        coo_data = arrow_table.column("soma_data").to_numpy()
        coo_coords = np.array(
            [
                arrow_table.column(f"soma_dim_{n}").to_numpy()
                for n in range(len(self.shape))
            ]
        ).T
        return pa.SparseCOOTensor.from_numpy(coo_data, coo_coords, shape=self.shape)


class SparseCSRMatrixReadIter(SparseTensorReadIterBase[pa.SparseCSRMatrix]):
    """Iterator over Arrow SparseCSRMatrix elements"""

    def __init__(self, sr: clib.SOMAReader, shape: NTuple):
        if len(shape) != 2:
            raise ValueError("CSR matrix format only supported for 2D SparseNDArray")
        super().__init__(sr, shape)

    def _from_table(self, arrow_table: pa.Table) -> pa.SparseCSRMatrix:
        if len(self.shape) != 2:
            raise SOMAError(f"internal error: expected shape == 2; got {self.shape}")
        data = arrow_table.column("soma_data").to_numpy()
        row = arrow_table.column("soma_dim_0").to_numpy()
        col = arrow_table.column("soma_dim_1").to_numpy()
        return pa.SparseCSRMatrix.from_scipy(
            sp.csr_matrix((data, (row, col)), shape=self.shape)
        )


class SparseCSCMatrixReadIter(SparseTensorReadIterBase[pa.SparseCSCMatrix]):
    """Iterator over Arrow SparseCSCMatrix elements"""

    def __init__(self, sr: clib.SOMAReader, shape: NTuple):
        if len(shape) != 2:
            raise ValueError("CSC matrix format only supported for 2D SparseNDArray")
        super().__init__(sr, shape)

    def _from_table(self, arrow_table: pa.Table) -> pa.SparseCSCMatrix:
        if len(self.shape) != 2:
            raise SOMAError(f"internal error: expected shape == 2; got {self.shape}")
        data = arrow_table.column("soma_data").to_numpy()
        row = arrow_table.column("soma_dim_0").to_numpy()
        col = arrow_table.column("soma_dim_1").to_numpy()
        return pa.SparseCSCMatrix.from_scipy(
            sp.csc_matrix((data, (row, col)), shape=self.shape)
        )
