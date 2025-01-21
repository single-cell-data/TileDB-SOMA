# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

from __future__ import annotations

import collections.abc
import math
from typing import Any, List, Literal, Sequence, Tuple, Union, cast

import numpy as np
import numpy.typing as npt
import pyarrow as pa
import scipy.sparse
from typing_extensions import TypeAlias

from .options._soma_tiledb_context import SOMATileDBContext
from .pytiledbsoma.fastercsx import compress_coo, copy_csx_to_dense, sort_csx_indices

NDArrayIndex: TypeAlias = npt.NDArray[np.integer[Any]]
NDArrayNumber: TypeAlias = npt.NDArray[Union[np.integer[Any], np.floating[Any]]]


class CompressedMatrix:
    """Fast compressed matrix construction from COO data.

    Private class intended _only_ for compressed matrix construction, for package-internal use.
    Export constructed matrix to SciPy (or comparable package) for use of the result.
    """

    __slots__ = (
        "indptr",
        "indices",
        "data",
        "shape",
        "format",
        "is_sorted",
        "no_duplicates",
        "context",
    )

    def __init__(
        self,
        indptr: NDArrayIndex,
        indices: NDArrayIndex,
        data: NDArrayNumber,
        shape: Tuple[int, int],
        format: Literal["csc", "csr"],
        is_sorted: bool,
        no_duplicates: bool | None,
        context: SOMATileDBContext,
    ) -> None:
        """Construct from PJV format. Not intended for direct use - use instead the
        static factory methods `from_ijd`, `from_pjv` and `from_soma`.
        """
        assert len(data) == len(indices)
        assert len(data) <= np.iinfo(indptr.dtype).max
        assert shape[1] <= np.iinfo(indices.dtype).max
        assert indptr[-1] == len(data) and indptr[0] == 0

        self.shape = shape
        self.indptr = indptr
        self.indices = indices
        self.data = data
        self.format = format
        self.is_sorted = is_sorted
        self.no_duplicates = no_duplicates
        self.context = context

    @staticmethod
    def from_ijd(
        i: NDArrayIndex | Sequence[NDArrayIndex],
        j: NDArrayIndex | Sequence[NDArrayIndex],
        d: NDArrayNumber | Sequence[NDArrayNumber],
        shape: Tuple[int, int],
        format: Literal["csc", "csr"],
        make_sorted: bool,
        context: SOMATileDBContext,
    ) -> CompressedMatrix:
        """Factory method accepting COO points stored in IJD vectors."""
        i = i if isinstance(i, collections.abc.Sequence) else (i,)
        j = j if isinstance(j, collections.abc.Sequence) else (j,)
        d = d if isinstance(d, collections.abc.Sequence) else (d,)

        if format == "csc":
            i, j = j, i
            n_major, n_minor = shape[1], shape[0]
        else:
            n_major, n_minor = shape[0], shape[1]

        nnz = sum(len(d_) for d_ in d)
        index_dtype = CompressedMatrix._smallest_index_dtype(nnz)
        indptr = np.zeros((n_major + 1), dtype=index_dtype)
        indices = np.empty((nnz,), dtype=index_dtype)
        data = np.empty((nnz,), dtype=d[0].dtype)
        compress_coo(
            context.native_context, (n_major, n_minor), i, j, d, indptr, indices, data
        )
        no_duplicates = None  # aka, unknown
        if make_sorted:
            no_duplicates = sort_csx_indices(
                context.native_context, indptr, indices, data
            )
        return CompressedMatrix(
            indptr, indices, data, shape, format, make_sorted, no_duplicates, context
        )

    @staticmethod
    def from_soma(
        tables: pa.Table | Sequence[pa.Table],
        shape: Tuple[int, int],
        format: Literal["csc", "csr"],
        make_sorted: bool,
        context: SOMATileDBContext,
    ) -> CompressedMatrix:
        """Factory method accepting a sequence of Arrow tables containing SOMA sparse matrix data.

        Table names must conform to the standard SOMA format, i.e., have three columns, named
        ``soma_dim_0``, ``soma_dim_1`` and ``soma_data``. All arrays in each table column must
        contain the same chunk count/size, and the dimension columns must be int64.
        """
        tbl = (
            pa.concat_tables(tables)
            if isinstance(tables, collections.abc.Sequence)
            else tables
        )

        def chunks(a: pa.Array | pa.ChunkedArray) -> List[pa.Array]:
            return (
                list(a) if isinstance(a, pa.Array) else cast(List[pa.Array], a.chunks)
            )

        if len(tbl) > 0:
            i = tuple(a.to_numpy() for a in chunks(tbl["soma_dim_0"]))
            j = tuple(a.to_numpy() for a in chunks(tbl["soma_dim_1"]))
            d = tuple(a.to_numpy() for a in chunks(tbl["soma_data"]))
        else:
            i = (tbl["soma_dim_0"].to_numpy(),)
            j = (tbl["soma_dim_1"].to_numpy(),)
            d = (tbl["soma_data"].to_numpy(),)

        return CompressedMatrix.from_ijd(i, j, d, shape, format, make_sorted, context)

    @property
    def nnz(self) -> int:
        return len(self.indices)

    @property
    def nbytes(self) -> int:
        return int(self.indptr.nbytes + self.indices.nbytes + self.data.nbytes)

    @property
    def dtype(self) -> npt.DTypeLike:
        return self.data.dtype

    def to_scipy(
        self, index: slice | None = None
    ) -> scipy.sparse.csr_matrix | scipy.sparse.csc_matrix:
        """Extract a slice on the compressed dimension and return as a
        :class:`scipy.sparse.csr_matrix` or
        :class:`scipy.sparse.csc_matrix`.

        Optionally allows slicing on compressed dimension during conversion, in which case
        an extra memory copy is avoided for the SOMA fast path (no-dups).
        """
        index = index or slice(None)
        assert isinstance(index, slice)
        assert index.step in (1, None)
        idx_start, idx_end, _ = index.indices(self.indptr.shape[0] - 1)
        n_major = max(idx_end - idx_start, 0)

        if n_major == self.indptr.shape[0] - 1:
            assert idx_start == 0 and n_major == idx_end
            indptr = self.indptr
            indices = self.indices
            data = self.data

        elif n_major == 0:
            return (
                scipy.sparse.csr_matrix((0, self.shape[1]), dtype=self.dtype)
                if self.format == "csr"
                else scipy.sparse.csc_matrix((self.shape[0], 0), dtype=self.dtype)
            )

        else:
            indptr = (self.indptr[idx_start : idx_end + 1]).copy()
            indices = self.indices[indptr[0] : indptr[-1]].copy()
            data = self.data[indptr[0] : indptr[-1]].copy()
            indptr -= indptr[0]

        shape = (
            (n_major, self.shape[1])
            if self.format == "csr"
            else (self.shape[0], n_major)
        )
        return CompressedMatrix._to_scipy(
            indptr,
            indices,
            data,
            shape,
            self.format,
            self.is_sorted,
            self.no_duplicates,
        )

    def to_numpy(self, index: slice | None = None) -> NDArrayNumber:
        """Extract a slice on the compressed dimension and return as a dense :class:`numpy.ndarray`."""
        index = index or slice(None)
        assert isinstance(index, slice)
        assert index.step in (1, None)
        major_idx_start, major_idx_end, _ = index.indices(self.indptr.shape[0] - 1)
        n_major = max(major_idx_end - major_idx_start, 0)

        out_shape = (
            (n_major, self.shape[1])
            if self.format == "csr"
            else (self.shape[0], n_major)
        )
        out = np.zeros(math.prod(out_shape), dtype=self.data.dtype)
        copy_csx_to_dense(
            self.context.native_context,
            major_idx_start,
            major_idx_end,
            self.shape,
            self.format,
            self.indptr,
            self.indices,
            self.data,
            out,
        )
        return out.reshape(out_shape)

    @classmethod
    def _smallest_index_dtype(cls, max_val: int) -> npt.DTypeLike:
        """NB: the underlying C++ code supports other index types, including uint16 and
        uint32. This helper method uses only int32/int64 to retain compatibility with
        the SciPy sparse matrix/array package.
        """
        candidate_index_types: list[npt.DTypeLike] = [np.int32, np.int64]
        for dt in candidate_index_types:
            if max_val <= np.iinfo(dt).max:
                return dt

        raise ValueError(
            "Unable to find index type sufficiently large for max index value."
        )

    @staticmethod
    def _to_scipy(
        indptr: NDArrayNumber,
        indices: NDArrayNumber,
        data: NDArrayNumber,
        shape: Tuple[int, int],
        format: Literal["csc", "csr"],
        is_sorted: bool,
        no_duplicates: bool | None,
    ) -> scipy.sparse.csr_matrix | scipy.sparse.csc_matrix:
        """
        This is to bypass the O(N) scan that :meth:`sparse._cs_matrix.__init__`
        performs when a new compressed matrix is created.

        See `SciPy bug 11496 <https://github.com/scipy/scipy/issues/11496>` for details.

        Conceptually, this is identical to:
            sparse.csr_matrix((data, indices, indptr), shape=shape)
        (or the equivalent for csc_matrix)
        """
        if is_sorted and bool(no_duplicates):
            # Fast path for SOMA common case.
            matrix = (
                scipy.sparse.csr_matrix.__new__(scipy.sparse.csr_matrix)
                if format == "csr"
                else scipy.sparse.csc_matrix.__new__(scipy.sparse.csc_matrix)
            )
            matrix.data = data
            matrix.indptr = indptr
            matrix.indices = indices
            matrix._shape = shape
            matrix._has_sorted_values = True
            matrix._has_canonical_format = True
            return matrix

        else:
            return (
                scipy.sparse.csr_matrix((data, indices, indptr), shape=shape)
                if format == "csr"
                else scipy.sparse.csc_matrix((data, indices, indptr), shape=shape)
            )
