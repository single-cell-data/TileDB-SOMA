import os
from concurrent.futures import Executor, ThreadPoolExecutor, wait
from typing import Callable, List, NamedTuple, Tuple, Type, Union, cast

import numba
import numba.typed
import numpy as np
import numpy.typing as npt
import pyarrow as pa
import scipy.sparse as sp
from somacore.query.types import IndexFactory, IndexLike

from ._eager_iter import EagerIterator
from ._funcs import _T, _Params
from ._sparse_nd_array import SparseNDArray
from ._types import NPIntArray, NPNDArray

try:
    # We need to `typeguard_ignore` the `@numba.jit`'d functions later in this file.
    # However, for some reason `from ._funcs import typeguard_ignore` does not work here: tests raise errors like:
    # ```
    # E numba.core.errors.TypingError: Failed in nopython mode pipeline (step: nopython frontend)
    # E Untyped global name 'TypeCheckMemo': Cannot determine Numba type of <class 'type'>
    # E
    # E File "src/tiledbsoma/_fast_csr.py", line 23:
    # E
    # E def read_csr(
    # E ^
    # ```
    # Directly importing here is the only known solution.
    from typeguard import typeguard_ignore
except ImportError:
    # Define a typeguard_ignore function so that we can use the `@typeguard_ignore`
    # decorator without having to depend upon typeguard at runtime.
    def typeguard_ignore(f: Callable[_Params, _T]) -> Callable[_Params, _T]:
        """No-op. Returns the argument unchanged."""
        return f


def read_csr(
    matrix: SparseNDArray,
    obs_joinids: pa.IntegerArray,
    var_joinids: pa.IntegerArray,
    index_factory: IndexFactory,
) -> "AccumulatedCSR":
    if not isinstance(matrix, SparseNDArray) or matrix.ndim != 2:
        raise TypeError("Can only read from a 2D SparseNDArray")

    max_workers = (os.cpu_count() or 4) + 2
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        acc = _CSRAccumulator(
            obs_joinids=obs_joinids,
            var_joinids=var_joinids,
            pool=pool,
            index_factory=index_factory,
        )
        for tbl in EagerIterator(
            matrix.read((obs_joinids, var_joinids)).tables(),
            pool=pool,
        ):
            acc.append(tbl["soma_dim_0"], tbl["soma_dim_1"], tbl["soma_data"])

        return acc.finalize()


class AccumulatedCSR(NamedTuple):
    """
    Private.

    Return type for the _CSRAccumulator.finalize method.
    Contains a sparse CSR's constituent elements.
    """

    data: NPNDArray
    indptr: NPIntArray
    indices: NPIntArray
    shape: Tuple[int, int]

    def to_scipy(self) -> sp.csr_matrix:
        """Create a Scipy ``sparse.csr_matrix`` from component elements.

        Conceptually, this is identical to::

            sparse.csr_matrix((data, indices, indptr), shape=shape)

        This ugliness is to bypass the O(N) scan that
        :meth:`sparse._cs_matrix.__init__`
        does when a new compressed matrix is created.

        See `SciPy bug 11496 <https://github.com/scipy/scipy/issues/11496>`
        for details.
        """
        matrix = sp.csr_matrix.__new__(sp.csr_matrix)
        matrix.data = self.data
        matrix.indptr = self.indptr
        matrix.indices = self.indices
        matrix._shape = self.shape
        return matrix


class _CSRAccumulator:
    """
    Fast accumulator of a CSR, based upon COO input.
    """

    def __init__(
        self,
        obs_joinids: pa.IntegerArray,
        var_joinids: pa.IntegerArray,
        pool: Executor,
        index_factory: IndexFactory,
    ):
        self.obs_joinids = obs_joinids
        self.var_joinids = var_joinids
        self.pool = pool

        self.shape: Tuple[int, int] = (len(self.obs_joinids), len(self.var_joinids))
        self.obs_indexer = index_factory(self.obs_joinids)
        self.var_indexer = index_factory(self.var_joinids)
        self.row_length: NPIntArray = np.zeros(
            (self.shape[0],), dtype=_select_dtype(self.shape[1])
        )

        # COO accumulated chunks, stored as list of triples (row_ind, col_ind, data)
        self.coo_chunks: List[
            Tuple[
                NPIntArray,  # row_ind
                NPIntArray,  # col_ind
                NPNDArray,  # data
            ]
        ] = []

    def append(
        self,
        row_joinids: Union[pa.Array, pa.ChunkedArray],
        col_joinids: Union[pa.Array, pa.ChunkedArray],
        data: Union[pa.Array, pa.ChunkedArray],
    ) -> None:
        """
        At accumulation time, do several things:

        * re-index to positional indices, and if possible, cast to smaller dtype
          to minimize memory footprint (at cost of some amount of time)
        * accumulate column counts by row, i.e., build the basis of the indptr
        * cache the tuple of data, row, col
        """
        rows_future = self.pool.submit(
            _reindex_and_cast,
            self.obs_indexer,
            row_joinids.to_numpy(),
            _select_dtype(self.shape[0]),
        )
        cols_future = self.pool.submit(
            _reindex_and_cast,
            self.var_indexer,
            col_joinids.to_numpy(),
            _select_dtype(self.shape[1]),
        )
        row_ind = rows_future.result()
        col_ind = cols_future.result()
        self.coo_chunks.append((row_ind, col_ind, data.to_numpy()))  # type: ignore[arg-type]
        _accum_row_length(self.row_length, row_ind)

    def finalize(self) -> AccumulatedCSR:
        nnz = sum(len(chunk[2]) for chunk in self.coo_chunks)
        index_dtype = _select_dtype(nnz)
        if nnz == 0:
            # There is no way to infer matrix dtype, so use a default and return
            # an empty matrix. Float32 is used as a default type, as it is most
            # compatible with AnnData expectations.
            empty = sp.csr_matrix((0, 0), dtype=np.float32)
            return AccumulatedCSR(
                data=empty.data,
                indptr=empty.indptr,
                indices=empty.indices,
                shape=self.shape,
            )

        # cumsum row lengths to get indptr
        indptr: NPIntArray = np.empty((self.shape[0] + 1,), dtype=index_dtype)
        indptr[0:1] = 0
        np.cumsum(self.row_length, out=indptr[1:])

        # Parallel copy of data and column indices
        indices: NPIntArray = np.empty((nnz,), dtype=index_dtype)
        data: NPNDArray = np.empty((nnz,), dtype=self.coo_chunks[0][2].dtype)

        # Empirically determined value. Needs to be large enough for reasonable
        # concurrency, without excessive write cache conflict. Controls the
        # number of rows that are processed in a single thread, and therefore
        # is the primary tuning parameter related to concurrency.
        row_rng_mask_bits = 18

        n_jobs = (self.shape[0] >> row_rng_mask_bits) + 1
        chunk_list = numba.typed.List(self.coo_chunks)
        wait(
            [
                self.pool.submit(
                    _copy_chunklist_range,
                    chunk_list,
                    data,
                    indices,
                    indptr,
                    row_rng_mask_bits,
                    job,
                )
                for job in range(n_jobs)
            ]
        )
        _finalize_indptr(indptr)
        return AccumulatedCSR(
            data=data, indptr=indptr, indices=indices, shape=self.shape
        )


@typeguard_ignore  # type: ignore[misc]
@numba.jit(nopython=True, nogil=True)  # type: ignore[attr-defined,misc]
def _accum_row_length(
    row_length: npt.NDArray[np.int64], row_ind: npt.NDArray[np.int64]
) -> None:
    for rind in row_ind:
        row_length[rind] += 1


@typeguard_ignore  # type: ignore[misc]
@numba.jit(nopython=True, nogil=True)  # type: ignore[attr-defined,misc]
def _copy_chunk_range(
    row_ind_chunk: npt.NDArray[np.signedinteger[npt.NBitBase]],
    col_ind_chunk: npt.NDArray[np.signedinteger[npt.NBitBase]],
    data_chunk: NPNDArray,
    data: NPNDArray,
    indices: npt.NDArray[np.signedinteger[npt.NBitBase]],
    indptr: npt.NDArray[np.signedinteger[npt.NBitBase]],
    row_rng_mask: int,
    row_rng_val: int,
) -> None:
    for n in range(len(data_chunk)):
        row = row_ind_chunk[n]
        if (row & row_rng_mask) != row_rng_val:
            continue
        ptr = indptr[row]
        indices[ptr] = col_ind_chunk[n]
        data[ptr] = data_chunk[n]
        indptr[row] += 1


@typeguard_ignore  # type: ignore[misc]
@numba.jit(nopython=True, nogil=True)  # type: ignore[attr-defined,misc]
def _copy_chunklist_range(
    chunk_list: numba.typed.List,
    data: NPNDArray,
    indices: npt.NDArray[np.signedinteger[npt.NBitBase]],
    indptr: npt.NDArray[np.signedinteger[npt.NBitBase]],
    row_rng_mask_bits: int,
    job: int,
) -> None:
    assert row_rng_mask_bits >= 1 and row_rng_mask_bits < 64
    row_rng_mask = (2**64 - 1) >> row_rng_mask_bits << row_rng_mask_bits
    row_rng_val = job << row_rng_mask_bits
    for row_ind_chunk, col_ind_chunk, data_chunk in chunk_list:
        _copy_chunk_range(
            row_ind_chunk,
            col_ind_chunk,
            data_chunk,
            data,
            indices,
            indptr,
            row_rng_mask,
            row_rng_val,
        )


@typeguard_ignore  # type: ignore[misc]
@numba.jit(nopython=True, nogil=True)  # type: ignore[attr-defined,misc]
def _finalize_indptr(indptr: npt.NDArray[np.signedinteger[npt.NBitBase]]) -> None:
    prev = 0
    for r in range(len(indptr)):
        t = indptr[r]
        indptr[r] = prev
        prev = t


def _select_dtype(
    maxval: int,
) -> Union[Type[np.int32], Type[np.int64]]:
    """
    Ascertain the "best" dtype for a zero-based index. Given our
    goal of minimizing memory use, "best" is currently defined as
    smallest.
    """
    if maxval > np.iinfo(np.int32).max:
        return np.int64
    else:
        return np.int32


def _reindex_and_cast(
    index: IndexLike, ids: npt.NDArray[np.int64], target_dtype: npt.DTypeLike
) -> npt.NDArray[np.int64]:
    return cast(
        npt.NDArray[np.int64], index.get_indexer(ids).astype(target_dtype, copy=False)
    )
