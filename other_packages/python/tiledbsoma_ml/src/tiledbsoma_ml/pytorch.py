# Copyright (c) 2021-2024 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2024 TileDB, Inc.
#
# Licensed under the MIT License.

from __future__ import annotations

import gc
import itertools
import logging
import os
import sys
import time
from contextlib import contextmanager
from itertools import islice
from math import ceil
from typing import (
    TYPE_CHECKING,
    Any,
    Dict,
    Iterable,
    Iterator,
    Sequence,
    Tuple,
    TypeVar,
    Union,
    cast,
)

import attrs
import numba
import numpy as np
import numpy.typing as npt
import pandas as pd
import pyarrow as pa
import scipy.sparse as sparse
import torch
import torchdata
from somacore.query._eager_iter import EagerIterator as _EagerIterator

import tiledbsoma as soma

logger = logging.getLogger("tiledbsoma_ml.pytorch")

_T = TypeVar("_T")
_T_co = TypeVar("_T_co", covariant=True)

if TYPE_CHECKING:
    NDArrayNumber = npt.NDArray[np.number[Any]]
else:
    NDArrayNumber = np.ndarray
XObsDatum = Tuple[NDArrayNumber, pd.DataFrame]
XObsNpDatum = Tuple[NDArrayNumber, NDArrayNumber]
XObsTensorDatum = Tuple[torch.Tensor, torch.Tensor]
"""Return type of ``ExperimentAxisQueryDataPipe`` that pairs a Tensor of ``obs`` row(s) with a Tensor of ``X`` matrix row(s).
The Tensors are rank 1 if ``batch_size`` is 1, otherwise the Tensors are rank 2."""


@attrs.define(frozen=True, kw_only=True)
class _ExperimentLocator:
    """State required to open the Experiment.

    Necessary as we will likely be invoked across multiple processes.

    Private.
    """

    uri: str
    tiledb_timestamp_ms: int
    tiledb_config: Dict[str, Union[str, float]]

    @classmethod
    def create(cls, experiment: soma.Experiment) -> "_ExperimentLocator":
        return _ExperimentLocator(
            uri=experiment.uri,
            tiledb_timestamp_ms=experiment.tiledb_timestamp_ms,
            tiledb_config=experiment.context.tiledb_config,
        )

    @contextmanager
    def open_experiment(self) -> Iterator[soma.Experiment]:
        context = soma.SOMATileDBContext(tiledb_config=self.tiledb_config)
        with soma.Experiment.open(
            self.uri, tiledb_timestamp=self.tiledb_timestamp_ms, context=context
        ) as exp:
            yield exp


class ExperimentAxisQueryIterable(Iterable[XObsDatum]):
    r"""Private base class for Dataset/DataPipe subclasses."""

    # XXX TODO - docstrings, slots, etc.

    def __init__(
        self,
        experiment: soma.Experiment,
        measurement_name: str = "RNA",
        X_name: str = "raw",
        obs_query: soma.AxisQuery | None = None,
        var_query: soma.AxisQuery | None = None,
        obs_column_names: Sequence[str] = ("soma_joinid",),
        batch_size: int = 1,
        shuffle: bool = True,
        io_batch_size: int = 2**17,
        shuffle_chunk_size: int = 64,
        seed: int | None = None,
        use_eager_fetch: bool = True,
        partition: bool = True,
    ):
        super().__init__()

        # Anything set in the instance needs to be picklable for multi-process users
        self.experiment_locator = _ExperimentLocator.create(experiment)
        self.measurement_name = measurement_name
        self.layer_name = X_name
        self.obs_query = obs_query
        self.var_query = var_query
        self.obs_column_names = list(obs_column_names)
        self.batch_size = batch_size
        self.io_batch_size = io_batch_size
        self.shuffle = shuffle
        self.use_eager_fetch = use_eager_fetch
        self._obs_joinids: npt.NDArray[np.int64] | None = None
        self._var_joinids: npt.NDArray[np.int64] | None = None
        self._shuffle_rng = np.random.default_rng(seed) if shuffle else None
        self.partition = partition
        self._initialized = False

        self.shuffle_chunk_size = shuffle_chunk_size
        if self.shuffle:
            # round io_batch_size up to a unit of shuffle_chunk_size to simplify code.
            self.io_batch_size = (
                ceil(io_batch_size / shuffle_chunk_size) * shuffle_chunk_size
            )

        if not self.obs_column_names:
            raise ValueError("Must specify at least one value in `obs_column_names`")

    def _create_obs_joinid_iter(self) -> Iterator[npt.NDArray[np.int64]]:
        """Private - create iterator over obs id chunks with split size of (roughly) io_batch_size.

        As appropriate, will chunk, shuffle and apply partitioning per worker.
        """
        assert self._obs_joinids is not None
        obs_joinids: npt.NDArray[np.int64] = self._obs_joinids

        if self.shuffle:
            assert self._shuffle_rng is not None
            assert self.io_batch_size % self.shuffle_chunk_size == 0
            shuffle_split = np.array_split(
                obs_joinids, max(1, ceil(len(obs_joinids) / self.shuffle_chunk_size))
            )
            self._shuffle_rng.shuffle(shuffle_split)
            obs_joinids_chunked = list(
                np.concatenate(b)
                for b in _batched(
                    shuffle_split, self.io_batch_size // self.shuffle_chunk_size
                )
            )
        else:
            obs_joinids_chunked = np.array_split(
                obs_joinids, max(1, ceil(len(obs_joinids) / self.io_batch_size))
            )

        # Now extract the partition for this worker
        partition, num_partitions = (
            _get_torch_partition_info() if self.partition else (0, 1)
        )
        obs_splits = _splits(len(obs_joinids_chunked), num_partitions)
        obs_partition_joinids = obs_joinids_chunked[
            obs_splits[partition] : obs_splits[partition + 1]
        ].copy()
        obs_joinid_iter = iter(obs_partition_joinids)

        if logger.isEnabledFor(logging.DEBUG) and self.partition:
            logger.debug(
                f"Process {os.getpid()} handling partition {partition + 1} of {num_partitions}, "
                f"partition_size={sum([len(chunk) for chunk in obs_partition_joinids])}"
            )

        return obs_joinid_iter

    def _init_once(self) -> None:
        """One-time per worker initialization. All operations be idempotent in order to support pipe reset()."""
        if self._initialized:
            return

        logger.debug(
            f"Initializing ExperimentAxisQueryIterable (shuffle={self.shuffle})"
        )

        with self.experiment_locator.open_experiment() as exp:
            query = exp.axis_query(
                measurement_name=self.measurement_name,
                obs_query=self.obs_query,
                var_query=self.var_query,
            )
            self._obs_joinids = query.obs_joinids().to_numpy()
            self._var_joinids = query.var_joinids().to_numpy()

        self._initialized = True

    def __iter__(self) -> Iterator[XObsDatum]:
        with self.experiment_locator.open_experiment() as exp:
            self._init_once()

            X = exp.ms[self.measurement_name].X[self.layer_name]
            if not isinstance(X, soma.SparseNDArray):
                raise NotImplementedError(
                    "ExperimentAxisQueryIterDataPipe only supported on X layers which are of type SparseNDArray"
                )

            obs_joinid_iter = self._create_obs_joinid_iter()
            _mini_batch_iter = self._mini_batch_iter(exp.obs, X, obs_joinid_iter)
            if self.use_eager_fetch:
                _mini_batch_iter = _EagerIterator(
                    _mini_batch_iter, pool=exp.context.threadpool
                )

            yield from _mini_batch_iter

    def __len__(self) -> int:
        self._init_once()
        assert self._obs_joinids is not None

        div, rem = divmod(len(self._obs_joinids), self.batch_size)
        return div + bool(rem)

    @property
    def shape(self) -> Tuple[int, int]:
        """Get the shape of the data that will be returned by this :class:`tiledbsoma.ml.pytorch.ExperimentDataPipe`.
        This is the number of obs (cell) and var (feature) counts in the returned data. If used in multiprocessing mode
        (i.e. :class:`torch.utils.data.DataLoader` instantiated with num_workers > 0), the obs (cell) count will reflect
        the size of the partition of the data assigned to the active process.

        Returns:
            A 2-tuple of ``int``s, for obs and var counts, respectively.

        Lifecycle:
            experimental
        """
        self._init_once()
        assert self._obs_joinids is not None
        assert self._var_joinids is not None
        return len(self._obs_joinids), len(self._var_joinids)

    def __getitem__(self, index: int) -> XObsDatum:
        raise NotImplementedError("Can only be iterated")

    def _io_batch_iter(
        self,
        obs: soma.DataFrame,
        X: soma.SparseNDArray,
        obs_joinid_iter: Iterator[npt.NDArray[np.int64]],
    ) -> Iterator[Tuple[sparse.csr_array, pd.DataFrame]]:
        """Iterate over IO batches, i.e., SOMA query/read, producing a tuple of
        (X: csr_array, obs: DataFrame).

        obs joinids read are controlled by the obs_joinid_iter. Iterator results will
        be reindexed and shuffled (if shuffling enabled).

        Private.
        """
        assert self._var_joinids is not None

        obs_column_names = (
            list(self.obs_column_names)
            if "soma_joinid" in self.obs_column_names
            else ["soma_joinid", *self.obs_column_names]
        )
        var_indexer = soma.IntIndexer(self._var_joinids, context=X.context)

        for obs_coords in obs_joinid_iter:
            st_time = time.perf_counter()
            obs_shuffled_coords = (
                obs_coords
                if self._shuffle_rng is None
                else self._shuffle_rng.permuted(obs_coords)
            )
            logger.debug(
                f"Retrieving next SOMA IO batch of length {len(obs_coords)}..."
            )

            obs_io_batch = cast(
                pd.DataFrame,
                obs.read(coords=(obs_coords,), column_names=obs_column_names)
                .concat()
                .to_pandas()
                .set_index("soma_joinid")
                .reindex(obs_shuffled_coords, copy=False)
                .reset_index(),
            )
            obs_io_batch = obs_io_batch[self.obs_column_names]

            X_tbl_iter = X.read(coords=(obs_coords, self._var_joinids)).tables()
            X_tbl = pa.concat_tables(X_tbl_iter)
            obs_indexer = soma.IntIndexer(obs_shuffled_coords, context=X.context)

            i = obs_indexer.get_indexer(X_tbl["soma_dim_0"].to_numpy())
            j = var_indexer.get_indexer(X_tbl["soma_dim_1"].to_numpy())
            d = X_tbl["soma_data"].to_numpy()
            if len(i) < np.iinfo(np.int32).max:
                # downcast where able, as this saves considerable memory
                i = i.astype(np.int32)
                j = j.astype(np.int32)
            X_io_batch = sparse.csr_array(
                (d, (i, j)),
                shape=(len(obs_coords), len(self._var_joinids)),
                copy=False,
            )

            del X_tbl, X_tbl_iter, i, j, d
            del obs_coords, obs_shuffled_coords, obs_indexer
            gc.collect()
            tm = time.perf_counter() - st_time
            logger.debug(
                f"Retrieved SOMA IO batch, took {tm:.2f}sec, {X_io_batch.shape[0]/tm:0.1f} samples/sec"
            )
            yield X_io_batch, obs_io_batch

    def _mini_batch_iter(
        self,
        obs: soma.DataFrame,
        X: soma.SparseNDArray,
        obs_joinid_iter: Iterator[npt.NDArray[np.int64]],
    ) -> Iterator[XObsDatum]:
        """Break IO batches into shuffled mini-batch-sized chunks.

        Private.
        """
        assert self._obs_joinids is not None
        assert self._var_joinids is not None

        io_batch_iter = self._io_batch_iter(obs, X, obs_joinid_iter)
        if self.use_eager_fetch:
            io_batch_iter = _EagerIterator(io_batch_iter, pool=X.context.threadpool)

        mini_batch_size = self.batch_size
        result: Tuple[sparse.csr_array | sparse.csr_matrix, pd.DataFrame] | None = (
            None  # partial result
        )
        for X_io_batch, obs_io_batch in io_batch_iter:
            assert X_io_batch.shape[0] == obs_io_batch.shape[0]
            assert X_io_batch.shape[1] == len(self._var_joinids)
            iob_idx = 0  # current offset into io batch
            iob_len = X_io_batch.shape[0]

            while iob_idx < iob_len:
                if result is None:
                    # perform zero copy slice where possible
                    result = (
                        X_io_batch[iob_idx : iob_idx + mini_batch_size],
                        obs_io_batch.iloc[iob_idx : iob_idx + mini_batch_size],
                    )
                    iob_idx += len(result[1])
                else:
                    # use remanent from previous IO batch
                    to_take = min(mini_batch_size - len(result[1]), iob_len - iob_idx)
                    result = (
                        # In older versions of scipy.sparse, vstack will return _matrix when
                        # called with _array. Various code paths must accomadate this bug (mostly
                        # in their type declarations)
                        sparse.vstack([result[0], X_io_batch[0:to_take]]),
                        pd.concat([result[1], obs_io_batch.iloc[0:to_take]]),
                    )
                    iob_idx += to_take

                assert result[0].shape[0] == result[1].shape[0]
                if result[0].shape[0] == mini_batch_size:
                    # yield result
                    yield (_csr_to_dense(result[0]), result[1])
                    result = None

        else:
            # yield a remnant, if any
            if result is not None:
                # yield result
                yield (_csr_to_dense(result[0]), result[1])


class ExperimentAxisQueryDataPipe(
    torchdata.datapipes.iter.IterDataPipe[  # type:ignore[misc]
        torch.utils.data.dataset.Dataset[XObsDatum]
    ],
):
    # TODO: XXX docstrings

    def __init__(
        self,
        experiment: soma.Experiment,
        measurement_name: str = "RNA",
        X_name: str = "raw",
        obs_query: soma.AxisQuery | None = None,
        var_query: soma.AxisQuery | None = None,
        obs_column_names: Sequence[str] = ("soma_joinid",),
        batch_size: int = 1,
        shuffle: bool = True,
        seed: int | None = None,
        io_batch_size: int = 2**17,
        shuffle_chunk_size: int = 64,
        use_eager_fetch: bool = True,
        # encoders: List[Encoder] | None = None,
    ):
        super().__init__()
        self._exp_iter = ExperimentAxisQueryIterable(
            experiment=experiment,
            measurement_name=measurement_name,
            X_name=X_name,
            obs_query=obs_query,
            var_query=var_query,
            obs_column_names=obs_column_names,
            batch_size=batch_size,
            shuffle=shuffle,
            seed=seed,
            io_batch_size=io_batch_size,
            use_eager_fetch=use_eager_fetch,
            shuffle_chunk_size=shuffle_chunk_size,
            # ---
            partition=True,
        )

    def __iter__(self) -> Iterator[XObsDatum]:
        batch_size = self._exp_iter.batch_size
        for X, obs in self._exp_iter:
            if batch_size == 1:
                X = X[0]
            yield X, obs

    def __len__(self) -> int:
        return self._exp_iter.__len__()

    @property
    def shape(self) -> Tuple[int, int]:
        return self._exp_iter.shape


class ExperimentAxisQueryIterableDataset(
    torch.utils.data.IterableDataset[XObsDatum]  # type:ignore[misc]
):
    # TODO: XXX docstrings

    def __init__(
        self,
        experiment: soma.Experiment,
        measurement_name: str = "RNA",
        X_name: str = "raw",
        obs_query: soma.AxisQuery | None = None,
        var_query: soma.AxisQuery | None = None,
        obs_column_names: Sequence[str] = ("soma_joinid",),
        batch_size: int = 1,  # XXX add docstring noting values >1 will not work with default collator
        shuffle: bool = True,
        seed: int | None = None,
        io_batch_size: int = 2**17,
        shuffle_chunk_size: int = 64,
        use_eager_fetch: bool = True,
    ):
        super().__init__()
        self._exp_iter = ExperimentAxisQueryIterable(
            experiment=experiment,
            measurement_name=measurement_name,
            X_name=X_name,
            obs_query=obs_query,
            var_query=var_query,
            obs_column_names=obs_column_names,
            batch_size=batch_size,
            shuffle=shuffle,
            seed=seed,
            io_batch_size=io_batch_size,
            use_eager_fetch=use_eager_fetch,
            shuffle_chunk_size=shuffle_chunk_size,
            # ---
            partition=True,
        )

    def __iter__(self) -> Iterator[XObsDatum]:
        batch_size = self._exp_iter.batch_size
        for X, obs in self._exp_iter:
            if batch_size == 1:
                X = X[0]
            yield X, obs

    def __len__(self) -> int:
        return self._exp_iter.__len__()

    @property
    def shape(self) -> Tuple[int, int]:
        return self._exp_iter.shape


def _collate_noop(datum: _T) -> _T:
    """Disable collation in dataloader instance."""
    return datum


def experiment_dataloader(
    ds: torchdata.datapipes.iter.IterDataPipe | torch.utils.data.IterableDataset,
    num_workers: int = 0,
    **dataloader_kwargs: Any,
) -> torch.utils.data.DataLoader:
    # TODO: XXX docstrings

    # XXX why do we disallow collate_fn?

    unsupported_dataloader_args = [
        "shuffle",
        "batch_size",
        "sampler",
        "batch_sampler",
        "collate_fn",
    ]
    if set(unsupported_dataloader_args).intersection(dataloader_kwargs.keys()):
        raise ValueError(
            f"The {','.join(unsupported_dataloader_args)} DataLoader params are not supported"
        )

    if num_workers > 0:
        _init_multiprocessing()

    return torch.utils.data.DataLoader(
        ds,
        batch_size=None,  # batching is handled by upstream iterator
        num_workers=num_workers,
        collate_fn=_collate_noop,
        # shuffling is handled by upstream iterator
        shuffle=False,
        **dataloader_kwargs,
    )


def _splits(total_length: int, sections: int) -> npt.NDArray[np.intp]:
    """For `total_length` points, compute start/stop offsets that split the length into roughly equal sizes.

    A total_length of L, split into N sections, will return L%N sections of size L//N+1,
    and the remainder as size L//N.

    This results in is the same split as numpy.array_split, for an array of length L.

    Examples
    --------
    >>> _splits(10, 3)
    array([0,  4,  7, 10])
    >>> _splits(4, 2)
    array([0, 2, 4])
    """
    if sections <= 0:
        raise ValueError("number of sections must greater than 0.") from None
    each_section, extras = divmod(total_length, sections)
    per_section_sizes = (
        [0] + extras * [each_section + 1] + (sections - extras) * [each_section]
    )
    splits = np.array(per_section_sizes, dtype=np.intp).cumsum()
    return splits


if sys.version_info >= (3, 12):
    _batched = itertools.batched

else:

    def _batched(
        iterable: Iterable[_T_co], n: int, *, strict: bool = False
    ) -> Iterator[Tuple[_T_co, ...]]:
        """Same as the Python 3.12+ itertools.batched -- polyfill for old Python versions."""
        if n < 1:
            raise ValueError("n must be at least one")
        it = iter(iterable)
        while batch := tuple(islice(it, n)):
            if strict and len(batch) != n:
                raise ValueError("batched(): incomplete batch")
            yield batch


def _get_torch_partition_info() -> Tuple[int, int]:
    """Return this workers partition and total partition count as a tuple.

    Private. Used to partition the iterator in some cases.

    Examples
    --------
    >>> _get_torch_partition_info()
    (0, 1)

    """
    worker_info = torch.utils.data.get_worker_info()
    if worker_info is None:
        loader_partition, loader_partitions = 0, 1
    else:
        loader_partition = worker_info.id
        loader_partitions = worker_info.num_workers

    if not torch.distributed.is_initialized():
        dist_partition, num_dist_partitions = 0, 1
    else:
        dist_partition = torch.distributed.get_rank()
        num_dist_partitions = torch.distributed.get_world_size()

    total_partitions = num_dist_partitions * loader_partitions
    partition = dist_partition * loader_partitions + loader_partition

    return partition, total_partitions


def _init_multiprocessing() -> None:
    """Ensures use of "spawn" for starting child processes with multiprocessing.

    Forked processes are known to be problematic:
      https://pytorch.org/docs/stable/notes/multiprocessing.html#avoiding-and-fighting-deadlocks
    Also, CUDA does not support forked child processes:
      https://pytorch.org/docs/stable/notes/multiprocessing.html#cuda-in-multiprocessing

    """
    orig_start_method = torch.multiprocessing.get_start_method()
    if orig_start_method != "spawn":
        if orig_start_method:
            logger.warning(
                "switching torch multiprocessing start method from "
                f'"{torch.multiprocessing.get_start_method()}" to "spawn"'
            )
        torch.multiprocessing.set_start_method("spawn", force=True)


@numba.njit(nogil=True, parallel=True)
def _csr_to_dense_inner(indptr, indices, data, out):  # type:ignore[no-untyped-def]
    n_rows = out.shape[0]
    for i in numba.prange(n_rows):
        for j in range(indptr[i], indptr[i + 1]):
            out[i, indices[j]] = data[j]

    return out


def _csr_to_dense(sp: sparse.csr_array | sparse.csr_matrix) -> NDArrayNumber:
    """Fast, parallel, variant of scipy.sparse.csr_array.todense.

    Typically 4-8X faster, dending on host and size of array/matrix.
    """
    assert isinstance(sp, (sparse.csr_array, sparse.csr_matrix))
    return cast(
        NDArrayNumber,
        _csr_to_dense_inner(
            sp.indptr, sp.indices, sp.data, np.zeros(sp.shape, dtype=sp.dtype)
        ),
    )
