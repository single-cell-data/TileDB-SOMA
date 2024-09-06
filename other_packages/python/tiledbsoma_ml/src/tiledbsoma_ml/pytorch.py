# Copyright (c) 2021-2024 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2024 TileDB, Inc.
#
# Licensed under the MIT License.

from __future__ import annotations

import contextlib
import gc
import itertools
import logging
import math
import os
import sys
import time
from contextlib import contextmanager
from itertools import islice
from math import ceil
from typing import (
    TYPE_CHECKING,
    Any,
    ContextManager,
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
import torch
import torchdata
from somacore.query._eager_iter import EagerIterator as _EagerIterator
from typing_extensions import Self, TypeAlias

import tiledbsoma as soma

logger = logging.getLogger("tiledbsoma_ml.pytorch")

_T = TypeVar("_T")
_T_co = TypeVar("_T_co", covariant=True)

if TYPE_CHECKING:
    # Python 3.8 does not support subscripting types, so work-around by
    # restricting this to when we are running a type checker.  TODO: remove
    # the conditional when Python 3.8 support is dropped.
    NDArrayNumber: TypeAlias = npt.NDArray[np.number[Any]]
else:
    NDArrayNumber: TypeAlias = np.ndarray

XObsDatum: TypeAlias = Tuple[NDArrayNumber, pd.DataFrame]
"""Return type of ``ExperimentAxisQueryIterableDataset`` and ``ExperimentAxisQueryIterDataPipe``,
which pairs a :class:`numpy.ndarray` of ``X`` row(s) with a :class:`pandas.DataFrame` of ``obs`` row(s). If the
``batch_size`` is 1, the objects are of rank 1, else they are of rank 2."""


@attrs.define(frozen=True, kw_only=True)
class _ExperimentLocator:
    """State required to open the Experiment.

    Necessary as we will likely be invoked across multiple processes.

    Private implementation class.
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
    """An :class:`Iterator` which reads ``X`` and ``obs`` data from a :class:`tiledbsoma.Experiment`, as
    selected by a user-specified :class:`tiledbsoma.ExperimentAxisQuery`. Each step of the iterator
    produces equal sized ``X`` and ``obs`` data, in the form of a :class:``numpy.ndarray`` and
    :class:`pandas.DataFrame`.

    Private base class for subclasses of :class:`torch.utils.data.IterableDataset` and
    :class:`torchdata.datapipes.iter.IterDataPipe`. Refer to :class:`ExperimentAxisQueryIterableDataset`
    and `ExperimentAxisQueryDataPipe` for more details on usage.

    Lifecycle:
        experimental
    """

    def __init__(
        self,
        query: soma.ExperimentAxisQuery,
        X_name: str,
        obs_column_names: Sequence[str] = ("soma_joinid",),
        batch_size: int = 1,
        shuffle: bool = True,
        io_batch_size: int = 2**16,
        shuffle_chunk_size: int = 64,
        seed: int | None = None,
        use_eager_fetch: bool = True,
    ):
        """
        Construct a new ``ExperimentAxisQueryIterable``, suitable for use with :class:`torch.utils.data.DataLoader`.

        The resulting iterator will produce a 2-tuple containing associated slices of ``X`` and ``obs`` data, as
        a NumPy ``ndarray`` and a Pandas ``DataFrame`` respectively.

        Args:
            query:
                A :class:`tiledbsoma.ExperimentAxisQuery`, defining the data which will be iterated over.
            X_name:
                The name of the X layer to read.
            obs_column_names:
                The names of the ``obs`` columns to return. At least one column name must be specified.
                Default is ``('soma_joinid',)``.
            batch_size:
                The number of rows of ``X`` and ``obs`` data to return in each iteration. Defaults to ``1``. A value of
                ``1`` will result in :class:`torch.Tensor` of rank 1 being returns (a single row); larger values will
                result in :class:`torch.Tensor`\ s of rank 2 (multiple rows).

                Note that a ``batch_size`` of 1 allows this ``IterableDataset`` to be used with :class:`torch.utils.data.DataLoader`
                batching, but you will achieve higher performance by performing batching in this class, and setting the ``DataLoader``
                batch_size parameter to ``None``.
            shuffle:
                Whether to shuffle the ``obs`` and ``X`` data being returned. Defaults to ``True``.
            io_batch_size:
                The number of ``obs``/``X`` rows to retrieve when reading data from SOMA. This impacts two aspects of
                this class's behavior: 1) The maximum memory utilization, with larger values providing
                better read performance, but also requiring more memory; 2) The number of rows read prior to shuffling
                (see ``shuffle`` parameter for details). The default value of 131,072 provides high performance, but
                may need to be reduced in memory limited hosts (or where a large number of :class:`DataLoader` workers
                are employed).
            shuffle_chunk_size:
                The number of contiguous rows sampled, prior to concatenation and shuffling.
                Larger numbers correspond to more randomness per training batch.
                If ``shuffle == False``, this parameter is ignored.
            seed:
                The random seed used for shuffling. Defaults to ``None`` (no seed). This arguiment *must* be specified when using
                :class:`torch.nn.parallel.DistributedDataParallel` to ensure data partitions are disjoint across worker
                processes.
            use_eager_fetch:
                Fetch the next SOMA chunk of ``obs`` and ``X`` data immediately after a previously fetched SOMA chunk is made
                available for processing via the iterator. This allows network (or filesystem) requests to be made in
                parallel with client-side processing of the SOMA data, potentially improving overall performance at the
                cost of doubling memory utilization. Defaults to ``True``.

        Returns:
            An ``iterable``, which can be iterated over using the Python ``iter()`` statement, or passed directly to
            a :class:`torch.data.utils.DataLoader` instance.

        Raises:
            ``ValueError`` on various unsupported or malformed parameter values.

        Lifecycle:
            experimental
        """

        super().__init__()

        # Anything set in the instance needs to be picklable for multi-process DataLoaders
        self.experiment_locator = _ExperimentLocator.create(query.experiment)
        self.layer_name = X_name
        self.measurement_name = query.measurement_name
        self.obs_query = query._matrix_axis_query.obs
        self.var_query = query._matrix_axis_query.var
        self.obs_column_names = list(obs_column_names)
        self.batch_size = batch_size
        self.io_batch_size = io_batch_size
        self.shuffle = shuffle
        self.use_eager_fetch = use_eager_fetch
        self._obs_joinids: npt.NDArray[np.int64] | None = None
        self._var_joinids: npt.NDArray[np.int64] | None = None
        self._shuffle_rng = np.random.default_rng(seed) if shuffle else None
        self.shuffle_chunk_size = shuffle_chunk_size
        self._initialized = False

        if self.shuffle:
            # round io_batch_size up to a unit of shuffle_chunk_size to simplify code.
            self.io_batch_size = (
                ceil(io_batch_size / shuffle_chunk_size) * shuffle_chunk_size
            )

        if not self.obs_column_names:
            raise ValueError("Must specify at least one value in `obs_column_names`")

    def _create_obs_joinid_iter(self) -> Iterator[npt.NDArray[np.int64]]:
        """Create iterator over obs id chunks with split size of (roughly) io_batch_size.

        As appropriate, will chunk, shuffle and apply partitioning per worker.

        IMPORTANT: in any scenario using torch.distributed, where WORLD_SIZE > 1, this will
        always partition such that each process has the same number of samples. Where
        the number of obs_joinids is not evenly divisible by the number of processes,
        the number of joinids will be dropped (dropped ids can never exceed WORLD_SIZE-1).

        Abstractly, the steps taken:
        1. Split the joinids into WORLD_SIZE sections (aka number of GPUS in DDP)
        2. Trim the splits to be of equal length
        3. Chunk and optionally shuffle the chunks
        4. Partition by number of data loader workers (to not generate redundant batches
           in cases where the DataLoader is running with `n_workers>1`).

        Private method.
        """
        assert self._obs_joinids is not None
        obs_joinids: npt.NDArray[np.int64] = self._obs_joinids

        # 1. Get the split for the model replica/GPU
        world_size, rank = _get_distributed_world_rank()
        _gpu_splits = _splits(len(obs_joinids), world_size)
        _gpu_split = obs_joinids[_gpu_splits[rank] : _gpu_splits[rank + 1]]

        # 2. Trip to be all of equal length
        min_len = np.diff(_gpu_splits).min()
        assert 0 <= (np.diff(_gpu_splits).min() - min_len) <= 1
        _gpu_split = _gpu_split[:min_len]

        # 3. Chunk and optionally shuffle chunks
        if self.shuffle:
            assert self._shuffle_rng is not None
            assert self.io_batch_size % self.shuffle_chunk_size == 0
            shuffle_split = np.array_split(
                _gpu_split, max(1, ceil(len(_gpu_split) / self.shuffle_chunk_size))
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
                _gpu_split, max(1, ceil(len(_gpu_split) / self.io_batch_size))
            )

        # 4. Partition by DataLoader worker
        n_workers, worker_id = _get_worker_world_rank()
        obs_splits = _splits(len(obs_joinids_chunked), n_workers)
        obs_partition_joinids = obs_joinids_chunked[
            obs_splits[worker_id] : obs_splits[worker_id + 1]
        ].copy()
        obs_joinid_iter = iter(obs_partition_joinids)

        if logger.isEnabledFor(logging.DEBUG):
            logger.debug(
                f"Process {os.getpid()} rank={rank}, world_size={world_size}, worker_id={worker_id}, n_workers={n_workers}, "
                f"partition_size={sum([len(chunk) for chunk in obs_partition_joinids])}"
            )

        return obs_joinid_iter

    def _init_once(self, exp: soma.Experiment | None = None) -> None:
        """One-time per worker initialization.

        All operations be idempotent in order to support pipe reset().

        Private method.
        """
        if self._initialized:
            return

        logger.debug(
            f"Initializing ExperimentAxisQueryIterable (shuffle={self.shuffle})"
        )

        if exp is None:
            # If no user-provided Experiment, open/close it ourselves
            exp_cm: ContextManager[soma.Experiment] = (
                self.experiment_locator.open_experiment()
            )
        else:
            # else, it is caller responsibility to open/close the experiment
            exp_cm = contextlib.nullcontext(exp)

        with exp_cm as exp:
            with exp.axis_query(
                measurement_name=self.measurement_name,
                obs_query=self.obs_query,
                var_query=self.var_query,
            ) as query:
                self._obs_joinids = query.obs_joinids().to_numpy()
                self._var_joinids = query.var_joinids().to_numpy()

        self._initialized = True

    def __iter__(self) -> Iterator[XObsDatum]:
        """Create iterator over query.

        Returns:
            ``iterator``

        Lifecycle:
            experimental
        """
        with self.experiment_locator.open_experiment() as exp:
            self._init_once(exp)
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
        """Return approximate number of batches this iterable will produce. If run in the context of :class:`torch.distributed` or
        as a multi-process loader (i.e., :class:`torch.utils.data.DataLoader` instantiated with num_workers > 0), the obs (cell)
        count will reflect the size of the partition of the data assigned to the active process.

        See import caveats in the PyTorch
        [:class:`torch.utils.data.DataLoader`](https://pytorch.org/docs/stable/data.html#torch.utils.data.DataLoader)
        domentation regarding ``len(dataloader)``, which also apply to this class.

        Returns:
            An ``int``.

        Lifecycle:
            experimental
        """
        # self._init_once()
        # assert self._obs_joinids is not None
        # world_size, _ = _get_distributed_world_rank()
        # n_workers, _ = _get_worker_world_rank()
        # div, rem = divmod(len(self._obs_joinids) // world_size, self.batch_size)
        # return div + bool(rem)
        return self.shape[0]

    @property
    def shape(self) -> Tuple[int, int]:
        """Get the approximate shape of the data that will be returned by this :class:`tiledbsoma_ml.ExperimentAxisQueryIterable`.
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
        world_size, _ = _get_distributed_world_rank()
        n_workers, _ = _get_worker_world_rank()
        div, rem = divmod(
            len(self._obs_joinids) // world_size // n_workers, self.batch_size
        )
        return div + bool(rem), len(self._var_joinids)

    def __getitem__(self, index: int) -> XObsDatum:
        raise NotImplementedError(
            "``ExperimentAxisQueryIterable can only be iterated - does not support mapping"
        )

    def _io_batch_iter(
        self,
        obs: soma.DataFrame,
        X: soma.SparseNDArray,
        obs_joinid_iter: Iterator[npt.NDArray[np.int64]],
    ) -> Iterator[Tuple[_CSR_IO_Buffer, pd.DataFrame]]:
        """Iterate over IO batches, i.e., SOMA query/read, producing a tuple of
        (X: csr_array, obs: DataFrame).

        obs joinids read are controlled by the obs_joinid_iter. Iterator results will
        be reindexed and shuffled (if shuffling enabled).

        Private method.
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
            obs_indexer = soma.IntIndexer(obs_shuffled_coords, context=X.context)
            logger.debug(
                f"Retrieving next SOMA IO batch of length {len(obs_coords)}..."
            )

            # to maximize optty's for concurrency, when in eager_fetch mode,
            # create the X read iterator first, as the eager iterator will begin
            # the read-ahead immediately. Then proceed to fetch obs DataFrame.
            # This matters most on latent backing stores, e.g., S3.
            #
            X_tbl_iter: Iterator[pa.Table] = X.read(
                coords=(obs_coords, self._var_joinids)
            ).tables()

            def make_csr(
                X_tbl: pa.Table,
                obs_coords: npt.NDArray[np.int64],
                var_coords: npt.NDArray[np.int64],
                obs_indexer: soma.IntIndexer,
            ) -> _CSR_IO_Buffer:
                """This function provides a GC after we throw off (large) garbage."""
                m = _CSR_IO_Buffer.from_ijd(
                    obs_indexer.get_indexer(X_tbl["soma_dim_0"]),
                    var_indexer.get_indexer(X_tbl["soma_dim_1"]),
                    X_tbl["soma_data"].to_numpy(),
                    shape=(len(obs_coords), len(var_coords)),
                )
                gc.collect(generation=0)
                return m

            _csr_iter = (
                make_csr(X_tbl, obs_coords, self._var_joinids, obs_indexer)
                for X_tbl in X_tbl_iter
            )
            if self.use_eager_fetch:
                _csr_iter = _EagerIterator(_csr_iter, pool=X.context.threadpool)

            # Now that X read is potentially in progress (in eager mode), go fetch obs data
            #
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

            X_io_batch = _CSR_IO_Buffer.merge(tuple(_csr_iter))

            del obs_indexer, obs_coords, obs_shuffled_coords, _csr_iter
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

        Private method.
        """
        assert self._obs_joinids is not None
        assert self._var_joinids is not None

        io_batch_iter = self._io_batch_iter(obs, X, obs_joinid_iter)
        if self.use_eager_fetch:
            io_batch_iter = _EagerIterator(io_batch_iter, pool=X.context.threadpool)

        mini_batch_size = self.batch_size
        result: Tuple[NDArrayNumber, pd.DataFrame] | None = None
        for X_io_batch, obs_io_batch in io_batch_iter:
            assert X_io_batch.shape[0] == obs_io_batch.shape[0]
            assert X_io_batch.shape[1] == len(self._var_joinids)
            iob_idx = 0  # current offset into io batch
            iob_len = X_io_batch.shape[0]

            while iob_idx < iob_len:
                if result is None:
                    # perform zero copy slice where possible
                    result = (
                        X_io_batch.densified_slice(
                            slice(iob_idx, iob_idx + mini_batch_size)
                        ),
                        obs_io_batch.iloc[iob_idx : iob_idx + mini_batch_size],
                    )
                    iob_idx += len(result[1])
                else:
                    # use remanent from previous IO batch
                    to_take = min(mini_batch_size - len(result[1]), iob_len - iob_idx)
                    result = (
                        np.concatenate(
                            [result[0], X_io_batch.densified_slice(slice(0, to_take))]
                        ),
                        pd.concat([result[1], obs_io_batch.iloc[0:to_take]]),
                    )
                    iob_idx += to_take

                assert result[0].shape[0] == result[1].shape[0]
                if result[0].shape[0] == mini_batch_size:
                    yield result
                    result = None

        else:
            # yield a remnant, if any
            if result is not None:
                yield result


class ExperimentAxisQueryDataPipe(
    torchdata.datapipes.iter.IterDataPipe[  # type:ignore[misc]
        torch.utils.data.dataset.Dataset[XObsDatum]
    ],
):
    """A :class:`torch.utils.data.IterableDataset` implementation that loads from a :class:`tiledbsoma.SOMAExperiment`.

    This class is based upon the now-deprecated :class:`torchdata.datapipes` API, and should only be used for
    legacy code. See [GitHub issue #1196](https://github.com/pytorch/data/issues/1196) and the
    TorchData [README](https://github.com/pytorch/data/blob/v0.8.0/README.md) for more information.

    See :class:`tiledbsoma_ml.ExperimentAxisQueryIterableDataset` for more information on using this class.

    Lifecycle:
        deprecated
    """

    def __init__(
        self,
        query: soma.ExperimentAxisQuery,
        X_name: str = "raw",
        obs_column_names: Sequence[str] = ("soma_joinid",),
        batch_size: int = 1,
        shuffle: bool = True,
        seed: int | None = None,
        io_batch_size: int = 2**16,
        shuffle_chunk_size: int = 64,
        use_eager_fetch: bool = True,
    ):
        """
        See :class:`tiledbsoma_ml.ExperimentAxisQueryIterableDataset` for more information on using this class.

        Lifecycle:
            deprecated
        """
        super().__init__()
        self._exp_iter = ExperimentAxisQueryIterable(
            query=query,
            X_name=X_name,
            obs_column_names=obs_column_names,
            batch_size=batch_size,
            shuffle=shuffle,
            seed=seed,
            io_batch_size=io_batch_size,
            use_eager_fetch=use_eager_fetch,
            shuffle_chunk_size=shuffle_chunk_size,
        )

    def __iter__(self) -> Iterator[XObsDatum]:
        """
        See :class:`tiledbsoma_ml.ExperimentAxisQueryIterableDataset` for more information on using this class.

        Lifecycle:
            deprecated
        """
        batch_size = self._exp_iter.batch_size
        for X, obs in self._exp_iter:
            if batch_size == 1:
                X = X[0]
            yield X, obs

    def __len__(self) -> int:
        """
        See :class:`tiledbsoma_ml.ExperimentAxisQueryIterableDataset` for more information on using this class.

        Lifecycle:
            deprecated
        """
        return self._exp_iter.__len__()

    @property
    def shape(self) -> Tuple[int, int]:
        """
        See :class:`tiledbsoma_ml.ExperimentAxisQueryIterableDataset` for more information on using this class.

        Lifecycle:
            deprecated
        """
        return self._exp_iter.shape


class ExperimentAxisQueryIterableDataset(
    torch.utils.data.IterableDataset[XObsDatum]  # type:ignore[misc]
):
    """A :class:`torch.utils.data.IterableDataset` implementation that loads from a :class:`tiledbsoma.SOMAExperiment`.

    This class works seamlessly with :class:`torch.utils.data.DataLoader` to load ``obs`` and ``X`` data as
    specified by a SOMA :class:`tiledbsoma.ExperimentAxisQuery`, providing an iterator over batches of
    ``obs`` and ``X`` data. Each iteration will yield a tuple containing an NumPy ndarray and a Pandas DataFrame.

    For example:

    >>> import torch
    >>> import tiledbsoma
    >>> import tiledbsoma_ml
    >>> with tiledbsoma.Experiment.open("my_experiment_path") as exp:
            with exp.axis_query(measurement_name="RNA", obs_query=tiledbsoma.AxisQuery(value_filter="tissue_type=='lung'")) as query:
                ds = tiledbsoma_ml.ExperimentAxisQueryIterableDataset(ds)
                dataloader = torch.utils.data.DataLoader(ds)
    >>> data = next(iter(dataloader))
    >>> data
    (array([0., 0., 0., ..., 0., 0., 0.], dtype=float32),
    soma_joinid
    0     57905025)
    >>> data[0]
    array([0., 0., 0., ..., 0., 0., 0.], dtype=float32)
    >>> data[1]
    soma_joinid
    0     57905025

    The ``batch_size`` parameter controls the number of rows of ``obs`` and ``X`` data that are returned in each
    iteration. If the ``batch_size`` is 1, then each result will have rank 1, else it will have rank 2. A ``batch_size``
    of 1 is compatible with :class:`torch.utils.data.DataLoader`-implemented batching, but it will usually be more
    performant to create mini-batches using this class, and set the ``DataLoader`` batch size to `None`.

    The ``obs_column_names`` parameter determines the data columns that are returned in the ``obs`` DataFrame (the
    default is a single column, containing the ``soma_joinid`` for the ``obs`` dimension).

    The `io_batch_size` parameter determines the number of rows read, from which mini-batches are yielded. A
    larger value will increase total memory usage, and may reduce average read time per row.

    Shuffling support is enabled with the ``shuffle`` parameter, and will normally be more performant than using
    :class:`DataLoader` shuffling. The shuffling algorithm works as follows:

      1. Rows selected by the query are subdivided into groups of size ``shuffle_chunk_size``, aka a "shuffle chunk".
      2. A random selection of shuffle `chunks` is drawn and read as a single I/O buffer (of size ``io_buffer_size``).
      3. The entire I/O buffer is shuffled.

    Put another way, we read randomly selected groups of observations from across all query results, concatenate
    those into an I/O buffer, and shuffle the buffer before returning mini-batches. The randomness of the shuffle
    is therefore determiend by the ``io_buffer_size`` (number of rows read), and the ``shuffle_chunk_size``
    (number of rows in each draw).

    Lifecycle:
        experimental
    """

    def __init__(
        self,
        query: soma.ExperimentAxisQuery,
        X_name: str = "raw",
        obs_column_names: Sequence[str] = ("soma_joinid",),
        batch_size: int = 1,
        shuffle: bool = True,
        seed: int | None = None,
        io_batch_size: int = 2**16,
        shuffle_chunk_size: int = 64,
        use_eager_fetch: bool = True,
    ):
        """
        Construct a new ``ExperimentAxisQueryIterable``, suitable for use with :class:`torch.utils.data.DataLoader`.

        The resulting iterator will produce a 2-tuple containing associated slices of ``X`` and ``obs`` data, as
        a NumPy ``ndarray`` and a Pandas ``DataFrame`` respectively.

        Args:
            query:
                A :class:`tiledbsoma.ExperimentAxisQuery`, defining the data which will be iterated over.
            X_name:
                The name of the X layer to read.
            obs_column_names:
                The names of the ``obs`` columns to return. At least one column name must be specified.
                Default is ``('soma_joinid',)``.
            batch_size:
                The number of rows of ``X`` and ``obs`` data to return in each iteration. Defaults to ``1``. A value of
                ``1`` will result in :class:`torch.Tensor` of rank 1 being returns (a single row); larger values will
                result in :class:`torch.Tensor`\ s of rank 2 (multiple rows).

                Note that a ``batch_size`` of 1 allows this ``IterableDataset`` to be used with :class:`torch.utils.data.DataLoader`
                batching, but you will achieve higher performance by performing batching in this class, and setting the ``DataLoader``
                batch_size parameter to ``None``.
            shuffle:
                Whether to shuffle the ``obs`` and ``X`` data being returned. Defaults to ``True``.
            io_batch_size:
                The number of ``obs``/``X`` rows to retrieve when reading data from SOMA. This impacts two aspects of
                this class's behavior: 1) The maximum memory utilization, with larger values providing
                better read performance, but also requiring more memory; 2) The number of rows read prior to shuffling
                (see ``shuffle`` parameter for details). The default value of 131,072 provides high performance, but
                may need to be reduced in memory limited hosts (or where a large number of :class:`DataLoader` workers
                are employed).
            shuffle_chunk_size:
                The number of contiguous rows sampled, prior to concatenation and shuffling.
                Larger numbers correspond to more randomness per training batch.
                If ``shuffle == False``, this parameter is ignored.
            seed:
                The random seed used for shuffling. Defaults to ``None`` (no seed). This arguiment *must* be specified when using
                :class:`torch.nn.parallel.DistributedDataParallel` to ensure data partitions are disjoint across worker
                processes.
            use_eager_fetch:
                Fetch the next SOMA chunk of ``obs`` and ``X`` data immediately after a previously fetched SOMA chunk is made
                available for processing via the iterator. This allows network (or filesystem) requests to be made in
                parallel with client-side processing of the SOMA data, potentially improving overall performance at the
                cost of doubling memory utilization. Defaults to ``True``.

        Returns:
            An ``iterable``, which can be iterated over using the Python ``iter()`` statement, or passed directly to
            a :class:`torch.data.utils.DataLoader` instance.

        Raises:
            ``ValueError`` on various unsupported or malformed parameter values.

        Lifecycle:
            experimental

        """
        super().__init__()
        self._exp_iter = ExperimentAxisQueryIterable(
            query=query,
            X_name=X_name,
            obs_column_names=obs_column_names,
            batch_size=batch_size,
            shuffle=shuffle,
            seed=seed,
            io_batch_size=io_batch_size,
            use_eager_fetch=use_eager_fetch,
            shuffle_chunk_size=shuffle_chunk_size,
        )

    def __iter__(self) -> Iterator[XObsDatum]:
        """Create Iterator yielding tuples of :class:`numpy.ndarray` and :class:`pandas.DataFrame`.

        Returns:
            ``iterator``

        Lifecycle:
            experimental
        """
        batch_size = self._exp_iter.batch_size
        for X, obs in self._exp_iter:
            if batch_size == 1:
                X = X[0]
            yield X, obs

    def __len__(self) -> int:
        """Return approximate number of batches this iterable will produce.

        See import caveats in the PyTorch
        [:class:`torch.utils.data.DataLoader`](https://pytorch.org/docs/stable/data.html#torch.utils.data.DataLoader)
        domentation regarding ``len(dataloader)``, which also apply to this class.

        Returns:
            An ``int``.

        Lifecycle:
            experimental
        """
        return self._exp_iter.__len__()

    @property
    def shape(self) -> Tuple[int, int]:
        """Get the shape of the data that will be returned by this :class:`tiledbsoma_ml.ExperimentAxisQueryIterableDataset`.

        This is the number of obs (cell) and var (feature) counts in the returned data. If used in multiprocessing mode
        (i.e. :class:`torch.utils.data.DataLoader` instantiated with num_workers > 0), the obs (cell) count will reflect
        the size of the partition of the data assigned to the active process.

        Returns:
            A 2-tuple of ``int``s, for obs and var counts, respectively.

        Lifecycle:
            experimental
        """
        return self._exp_iter.shape


def experiment_dataloader(
    ds: torchdata.datapipes.iter.IterDataPipe | torch.utils.data.IterableDataset,
    # num_workers: int = 0,
    **dataloader_kwargs: Any,
) -> torch.utils.data.DataLoader:
    """Factory method for :class:`torch.utils.data.DataLoader`. This method can be used to safely instantiate a
    :class:`torch.utils.data.DataLoader` that works with :class:`tiledbsoma_ml.ExperimentAxisQueryIterableDataset`
    or :class:`tiledbsoma_ml.ExperimentAxisQueryIterDataPipe`.

    Several :class:`torch.utils.data.DataLoader` constructor parameters are not applicable, or are non-performant,
    when using loaders form this module, including ``shuffle``, ``batch_size``, ``sampler``, and ``batch_sampler``.
    Specifying any of these parameters will result in an error.

    Refer to ``https://pytorch.org/docs/stable/data.html#torch.utils.data.DataLoader`` for more information on
    :class:`torch.utils.data.DataLoader` parameters.

    Args:
        ds:
            A :class:`torch.utils.data.IterableDataset` or a :class:`torchdata.datapipes.iter.IterDataPipe`. May
            include chained data pipes.
        num_workers:
            How many subprocesses to use for data loading. 0 means that the data will be loaded in the main process. (default: 0)
        **dataloader_kwargs:
            Additional keyword arguments to pass to the :class:`torch.utils.data.DataLoader` constructor,
            except for ``shuffle``, ``batch_size``, ``sampler``, and ``batch_sampler``, which are not
            supported when data loaders in this module.

    Returns:
        A :class:`torch.utils.data.DataLoader`.

    Raises:
        ValueError: if any of the ``shuffle``, ``batch_size``, ``sampler``, or ``batch_sampler`` params
            are passed as keyword arguments.

    Lifecycle:
        experimental
    """
    unsupported_dataloader_args = [
        "shuffle",
        "batch_size",
        "sampler",
        "batch_sampler",
    ]
    if set(unsupported_dataloader_args).intersection(dataloader_kwargs.keys()):
        raise ValueError(
            f"The {','.join(unsupported_dataloader_args)} DataLoader parameters are not supported"
        )

    if dataloader_kwargs.get("num_workers", 0) > 0:
        _init_multiprocessing()

    if "collate_fn" not in dataloader_kwargs:
        dataloader_kwargs["collate_fn"] = _collate_noop

    return torch.utils.data.DataLoader(
        ds,
        batch_size=None,  # batching is handled by upstream iterator
        shuffle=False,  # shuffling is handled by upstream iterator
        **dataloader_kwargs,
    )


def _collate_noop(datum: _T) -> _T:
    """Noop collation for use with a dataloader instance.

    Private.
    """
    return datum


def _splits(total_length: int, sections: int) -> npt.NDArray[np.intp]:
    """For `total_length` points, compute start/stop offsets that split the length into roughly equal sizes.

    A total_length of L, split into N sections, will return L%N sections of size L//N+1,
    and the remainder as size L//N. This results in the same split as numpy.array_split,
    for an array of length L and sections N.

    Private.

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

    def _batched(iterable: Iterable[_T_co], n: int) -> Iterator[Tuple[_T_co, ...]]:
        """Same as the Python 3.12+ itertools.batched -- polyfill for old Python versions."""
        if n < 1:
            raise ValueError("n must be at least one")
        it = iter(iterable)
        while batch := tuple(islice(it, n)):
            yield batch


def _get_distributed_world_rank() -> Tuple[int, int]:
    """Return tuple containing equivalent of torch.distributed world size and rank."""
    if torch.distributed.is_initialized():
        world_size = torch.distributed.get_world_size()
        rank = torch.distributed.get_rank()
    else:
        # sometimes these are set even before torch.distributed is initialized, e.g., by torchrun
        world_size = int(os.environ.get("WORLD_SIZE", 1))
        rank = int(os.environ.get("RANK", 0))

    return world_size, rank


def _get_worker_world_rank() -> Tuple[int, int]:
    """Return number of dataloader workers and our worker rank/id"""
    worker_info = torch.utils.data.get_worker_info()
    if worker_info is None:
        return 1, 0
    return worker_info.num_workers, worker_info.id


def _init_multiprocessing() -> None:
    """Ensures use of "spawn" for starting child processes with multiprocessing.

    Forked processes are known to be problematic:
      https://pytorch.org/docs/stable/notes/multiprocessing.html#avoiding-and-fighting-deadlocks
    Also, CUDA does not support forked child processes:
      https://pytorch.org/docs/stable/notes/multiprocessing.html#cuda-in-multiprocessing

    Private.
    """
    orig_start_method = torch.multiprocessing.get_start_method()
    if orig_start_method != "spawn":
        if orig_start_method:
            logger.warning(
                "switching torch multiprocessing start method from "
                f'"{torch.multiprocessing.get_start_method()}" to "spawn"'
            )
        torch.multiprocessing.set_start_method("spawn", force=True)


class _CSR_IO_Buffer:
    """Implement a minimal CSR matrix with specific optimizations for use in this package.

    Operations supported are:
    * Incrementally build a CSR from COO, allowing overlapped I/O and CSR conversion for I/O batches,
      and a final "merge" step which combines the result.
    * Zero intermediate copy conversion of an arbitrary row slice to dense (ie., mini-batch extraction).
    * Parallel ops where it makes sense (construction, merge, etc)

    Overall is significantly faster, and uses less memory, than the equivalent scipy.sparse operations.
    """

    __slots__ = ("indptr", "indices", "data", "shape")

    def __init__(
        self,
        indptr: NDArrayNumber,
        indices: NDArrayNumber,
        data: NDArrayNumber,
        shape: Tuple[int, int],
    ) -> None:
        """Construct from PJV format."""
        assert len(data) == len(indices)
        assert len(data) <= np.iinfo(indptr.dtype).max
        assert shape[1] <= np.iinfo(indices.dtype).max
        assert indptr[-1] == len(data) and indptr[0] == 0

        self.shape = shape
        self.indptr = indptr
        self.indices = indices
        self.data = data

    @staticmethod
    def from_ijd(
        i: NDArrayNumber, j: NDArrayNumber, d: NDArrayNumber, shape: Tuple[int, int]
    ) -> _CSR_IO_Buffer:
        """Factory from COO"""
        nnz = len(d)
        indptr = np.zeros((shape[0] + 1), dtype=smallest_uint_dtype(nnz))
        indices = np.empty((nnz,), dtype=smallest_uint_dtype(shape[1]))
        data = np.empty((nnz,), dtype=d.dtype)
        _coo_to_csr_inner(shape[0], i, j, d, indptr, indices, data)
        return _CSR_IO_Buffer(indptr, indices, data, shape)

    @staticmethod
    def from_pjd(
        p: NDArrayNumber, j: NDArrayNumber, d: NDArrayNumber, shape: Tuple[int, int]
    ) -> _CSR_IO_Buffer:
        """Factory from CSR"""
        return _CSR_IO_Buffer(p, j, d, shape)

    @property
    def nnz(self) -> int:
        return len(self.indices)

    @property
    def nbytes(self) -> int:
        return int(self.indptr.nbytes + self.indices.nbytes + self.data.nbytes)

    @property
    def dtype(self) -> npt.DTypeLike:
        return self.data.dtype

    def densified_slice(self, row_index: slice) -> NDArrayNumber:
        assert isinstance(row_index, slice)
        assert row_index.step in (1, None)
        row_idx_start, row_idx_end, _ = row_index.indices(self.indptr.shape[0] - 1)
        n_rows = max(row_idx_end - row_idx_start, 0)
        out = np.zeros((n_rows, self.shape[1]), dtype=self.data.dtype)
        if n_rows >= 0:
            _csr_to_dense_inner(
                row_idx_start, n_rows, self.indptr, self.indices, self.data, out
            )
        return out

    @staticmethod
    def merge(mtxs: Sequence[_CSR_IO_Buffer]) -> _CSR_IO_Buffer:
        assert len(mtxs) > 0
        nnz = sum(m.nnz for m in mtxs)
        shape = mtxs[0].shape
        for m in mtxs[1:]:
            assert m.shape == mtxs[0].shape
            assert m.indices.dtype == mtxs[0].indices.dtype
        assert all(m.shape == shape for m in mtxs)

        indptr = np.sum(
            [m.indptr for m in mtxs], axis=0, dtype=smallest_uint_dtype(nnz)
        )
        indices = np.empty((nnz,), dtype=mtxs[0].indices.dtype)
        data = np.empty((nnz,), mtxs[0].data.dtype)

        _csr_merge_inner(
            tuple((m.indptr.astype(indptr.dtype), m.indices, m.data) for m in mtxs),
            indptr,
            indices,
            data,
        )
        return _CSR_IO_Buffer.from_pjd(indptr, indices, data, shape)

    def sort_indices(self) -> Self:
        """Sort indices, IN PLACE."""
        _csr_sort_indices(self.indptr, self.indices, self.data)
        return self


def smallest_uint_dtype(max_val: int) -> npt.DTypeLike:
    for dt in [np.uint16, np.uint32]:
        if max_val <= np.iinfo(dt).max:
            return dt
    else:
        return np.uint64


@numba.njit(nogil=True, parallel=True)  # type:ignore[misc]
def _csr_merge_inner(
    As: Tuple[Tuple[NDArrayNumber, NDArrayNumber, NDArrayNumber], ...],  # P,J,D
    Bp: NDArrayNumber,
    Bj: NDArrayNumber,
    Bd: NDArrayNumber,
) -> None:
    n_rows = len(Bp) - 1
    offsets = Bp.copy()
    for Ap, Aj, Ad in As:
        n_elmts = Ap[1:] - Ap[:-1]
        for n in numba.prange(n_rows):
            Bj[offsets[n] : offsets[n] + n_elmts[n]] = Aj[Ap[n] : Ap[n] + n_elmts[n]]
            Bd[offsets[n] : offsets[n] + n_elmts[n]] = Ad[Ap[n] : Ap[n] + n_elmts[n]]
        offsets[:-1] += n_elmts


@numba.njit(nogil=True, parallel=True)  # type:ignore[misc]
def _csr_to_dense_inner(
    row_idx_start: int,
    n_rows: int,
    indptr: NDArrayNumber,
    indices: NDArrayNumber,
    data: NDArrayNumber,
    out: NDArrayNumber,
) -> None:
    for i in numba.prange(row_idx_start, row_idx_start + n_rows):
        for j in range(indptr[i], indptr[i + 1]):
            out[i - row_idx_start, indices[j]] = data[j]


@numba.njit(nogil=True, parallel=True, inline="always")  # type:ignore[misc]
def _count_rows(n_rows: int, Ai: NDArrayNumber, Bp: NDArrayNumber) -> NDArrayNumber:
    """Private: parallel row count."""
    nnz = len(Ai)

    partition_size = 32 * 1024**2
    n_partitions = math.ceil(nnz / partition_size)
    if n_partitions > 1:
        counts = np.zeros((n_partitions, n_rows), dtype=Bp.dtype)
        for p in numba.prange(n_partitions):
            for n in range(p * partition_size, min(nnz, (p + 1) * partition_size)):
                row = Ai[n]
                counts[p, row] += 1

        Bp[:-1] = counts.sum(axis=0)
    else:
        for n in range(nnz):
            row = Ai[n]
            Bp[row] += 1

    return Bp


@numba.njit(nogil=True, parallel=True)  # type:ignore[misc]
def _coo_to_csr_inner(
    n_rows: int,
    Ai: NDArrayNumber,
    Aj: NDArrayNumber,
    Ad: NDArrayNumber,
    Bp: NDArrayNumber,
    Bj: NDArrayNumber,
    Bd: NDArrayNumber,
) -> None:
    nnz = len(Ai)

    _count_rows(n_rows, Ai, Bp)

    # cum sum to get the row index pointers (NOTE: starting with zero)
    cumsum = 0
    for n in range(n_rows):
        tmp = Bp[n]
        Bp[n] = cumsum
        cumsum += tmp
    Bp[n_rows] = nnz

    # Reorganize all of the data. Side-effect: pointers shifted (reversed in the
    # subsequent section).
    #
    # Method is concurrent (partioned by rows) if number of rows is greater
    # than 2**partition_bits. This partitioning scheme leverages the fact
    # that reads are much cheaper than writes.
    #
    # The code is equivalent to:
    #   for n in range(nnz):
    #       row = Ai[n]
    #       dst_row = Bp[row]
    #       Bj[dst_row] = Aj[n]
    #       Bd[dst_row] = Ad[n]
    #       Bp[row] += 1

    partition_bits = 13
    n_partitions = (n_rows + 2**partition_bits - 1) >> partition_bits
    for p in numba.prange(n_partitions):
        for n in range(nnz):
            row = Ai[n]
            if (row >> partition_bits) != p:
                continue
            dst_row = Bp[row]
            Bj[dst_row] = Aj[n]
            Bd[dst_row] = Ad[n]
            Bp[row] += 1

    # Shift the pointers by one slot (ie., start at zero)
    prev_ptr = 0
    for n in range(n_rows + 1):
        tmp = Bp[n]
        Bp[n] = prev_ptr
        prev_ptr = tmp


@numba.njit(nogil=True, parallel=True)  # type:ignore[misc]
def _csr_sort_indices(Bp: NDArrayNumber, Bj: NDArrayNumber, Bd: NDArrayNumber) -> None:
    """In-place sort of minor axis indices"""
    n_rows = len(Bp) - 1
    for r in numba.prange(n_rows):
        row_start = Bp[r]
        row_end = Bp[r + 1]
        order = np.argsort(Bj[row_start:row_end])
        Bj[row_start:row_end] = Bj[row_start:row_end][order]
        Bd[row_start:row_end] = Bd[row_start:row_end][order]
