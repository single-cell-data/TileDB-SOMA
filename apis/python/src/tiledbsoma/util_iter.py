import abc
from typing import TypeVar

import numpy as np
import pyarrow as pa
import scipy.sparse as sp

# This package's pybind11 code
import tiledbsoma.libtiledbsoma as clib

from . import somacore  # to be replaced by somacore package, when available
from .types import NTuple


class TableReadIter(somacore.ReadIter[pa.Table]):
    def __init__(self, sr: clib.SOMAReader):
        self.sr = sr

    def __next__(self) -> pa.Table:
        arrow_table = self.sr.read_next()
        if arrow_table is None:
            raise StopIteration

        return arrow_table

    def concat(self) -> pa.Table:
        return pa.concat_tables(self)

    def close(self) -> None:
        """TODO: Implement when there is a stateful handle that needs to be closed"""
        pass


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
        arrow_tables = pa.concat(self)
        return self._from_table(arrow_tables)

    def close(self) -> None:
        """TODO: Implement when there is a stateful handle that needs to be closed"""
        pass


class SparseCOOTensorReadIter(SparseTensorReadIterBase[pa.SparseCOOTensor]):
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
    def __init__(self, sr: clib.SOMAReader, shape: NTuple):
        if len(shape) != 2:
            raise ValueError("CSR matrix format only supported for 2D SparseNDArray")
        super().__init__(sr, shape)

    def _from_table(self, arrow_table: pa.Table) -> pa.SparseCSRMatrix:
        assert len(self.shape) == 2, "Internal error"
        data = arrow_table.column("soma_data").to_numpy()
        row = arrow_table.column("soma_dim_0").to_numpy()
        col = arrow_table.column("soma_dim_1").to_numpy()
        return pa.SparseCSRMatrix.from_scipy(
            sp.csr_matrix((data, (row, col)), shape=self.shape)
        )


class SparseCSCMatrixReadIter(SparseTensorReadIterBase[pa.SparseCSCMatrix]):
    def __init__(self, sr: clib.SOMAReader, shape: NTuple):
        if len(shape) != 2:
            raise ValueError("CSC matrix format only supported for 2D SparseNDArray")
        super().__init__(sr, shape)

    def _from_table(self, arrow_table: pa.Table) -> pa.SparseCSCMatrix:
        assert len(self.shape) == 2, "Internal error"
        data = arrow_table.column("soma_data").to_numpy()
        row = arrow_table.column("soma_dim_0").to_numpy()
        col = arrow_table.column("soma_dim_1").to_numpy()
        return pa.SparseCSCMatrix.from_scipy(
            sp.csc_matrix((data, (row, col)), shape=self.shape)
        )

    # def __init__(self, sr: clib.SOMAReader, shape: NTuple):
    #     self.sr = sr
    #     self.shape = shape

    # def __next__(self) -> pa.SparseCOOTensor:
    #     arrow_table = self.sr.read_next()
    #     if arrow_table is None:
    #         raise StopIteration

    #     """
    #     Temporary work-around for Arrow bug.

    #     In PyArrow 9.X and 10.X, there is a bug preventing the creation of "empty"
    #     (zero element) SparseCOOTensor objects.

    #     See https://issues.apache.org/jira/browse/ARROW-17933

    #     Just stop the iteration when we run out of results. The caller must be
    #     prepared to have a StopIteration, rather than an empty tensor, as the result of
    #     an empty query.

    #     TODO: remove this work-around when above bug is fixed in Arrow.
    #     """
    #     if arrow_table.num_rows == 0:
    #         raise StopIteration

    #     return _sparseCOOTensor_from_Table(arrow_table, self.shape)

    # def concat(self) -> pa.SparseCOOTensor:
    #     arrow_tables = []
    #     while True:
    #         tbl = self.sr.read_next()
    #         if tbl is None:
    #             break
    #         arrow_tables.append(tbl)

    #     # TODO: temp work-around. See note above about bug in Arrow SparseCOOTensor
    #     if len(arrow_tables) == 1 and arrow_tables[0].num_rows == 0:
    #         return None

    #     return _sparseCOOTensor_from_Table(pa.concat_tables(arrow_tables), self.shape)

    # def close(self) -> None:
    #     """TODO: Implement when there is a stateful handle that needs to be closed"""
    #     pass


# def _sparseCOOTensor_from_Table(
#     arrow_table: pa.Table, shape: NTuple
# ) -> pa.SparseCOOTensor:
#     coo_data = arrow_table.column("soma_data").to_numpy()
#     coo_coords = np.array(
#         [arrow_table.column(f"soma_dim_{n}").to_numpy() for n in range(len(shape))]
#     ).T
#     return pa.SparseCOOTensor.from_numpy(coo_data, coo_coords, shape=shape)


# def _sparseCSCMatrix_from_Table(
#     arrow_table: pa.Table, shape: NTuple
# ) -> pa.SparseCSCMatrix:
#     assert len(shape) == 2, "Internal error"
#     data = arrow_table.column("soma_data").to_numpy()
#     row = arrow_table.column("soma_dim_0").to_numpy()
#     col = arrow_table.column("soma_dim_1").to_numpy()
#     return pa.SparseCSCMatrix.from_scipy(sp.csc_matrix((data, (row, col)), shape=shape))


# def _sparseCSRMatrix_from_Table(
#     arrow_table: pa.Table, shape: NTuple
# ) -> pa.SparseCSRMatrix:
#     assert len(shape) == 2, "Internal error"
#     data = arrow_table.column("soma_data").to_numpy()
#     row = arrow_table.column("soma_dim_0").to_numpy()
#     col = arrow_table.column("soma_dim_1").to_numpy()
#     return pa.SparseCSRMatrix.from_scipy(sp.csr_matrix((data, (row, col)), shape=shape))
