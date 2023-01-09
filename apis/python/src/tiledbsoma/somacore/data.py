#
# Read types
#
import abc
from typing import Iterator, TypeVar

import pyarrow

_T = TypeVar("_T")


# Sparse reads are returned as an iterable structure:


class ReadIter(Iterator[_T], metaclass=abc.ABCMeta):
    """SparseRead result iterator allowing users to flatten the iteration."""

    # __iter__ is already implemented as `return self` in Iterator.
    # SOMA implementations must implement __next__.

    @abc.abstractmethod
    def concat(self) -> _T:
        """Returns all the requested data in a single operation.

        If some data has already been retrieved using ``next``, this will return
        the rest of the data after that is already returned.
        """
        raise NotImplementedError()


class SparseRead:
    """Intermediate type to choose result format when reading a sparse array.

    A query may not be able to return all of these formats. The concrete result
    may raise a ``NotImplementedError`` or may choose to raise a different
    exception (likely a ``TypeError``) containing more specific information
    about why the given format is not supported.
    """

    def coos(self) -> ReadIter[pyarrow.SparseCOOTensor]:
        raise NotImplementedError()

    def cscs(self) -> ReadIter[pyarrow.SparseCSCMatrix]:
        raise NotImplementedError()

    def csrs(self) -> ReadIter[pyarrow.SparseCSRMatrix]:
        raise NotImplementedError()

    def dense_tensors(self) -> ReadIter[pyarrow.Tensor]:
        raise NotImplementedError()

    def record_batches(self) -> ReadIter[pyarrow.RecordBatch]:
        raise NotImplementedError()

    def tables(self) -> ReadIter[pyarrow.Table]:
        raise NotImplementedError()
