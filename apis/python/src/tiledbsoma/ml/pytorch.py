# Copyright (c) 2021-2024 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2024 TileDB, Inc.
#
# Licensed under the MIT License.

from __future__ import annotations

import gc
import itertools
import logging
import os
import typing
from contextlib import contextmanager
from datetime import timedelta
from math import ceil
from time import time
from typing import Any, Dict, Iterator, List, Optional, Sequence, Tuple, Union

import numpy as np
import numpy.typing as npt
import pandas as pd
import psutil
import torch
import torchdata.datapipes.iter as pipes
from attrs import define
from numpy.random import Generator
from pyarrow import Table
from scipy import sparse
from somacore.query._eager_iter import EagerIterator as _EagerIterator
from torch import Tensor
from torch import distributed as dist
from torch.utils.data import DataLoader
from torch.utils.data.dataset import Dataset

import tiledbsoma as soma

from .encoders import Encoder, LabelEncoder

pytorch_logger = logging.getLogger("tiledbsoma.ml.pytorch")

# TODO: Rename to reflect the correct order of the Tensors within the tuple: (X, obs)
ObsAndXDatum = Tuple[Tensor, Tensor]
"""Return type of ``ExperimentDataPipe`` that pairs a Tensor of ``obs`` row(s) with a Tensor of ``X`` matrix row(s).
The Tensors are rank 1 if ``batch_size`` is 1, otherwise the Tensors are rank 2."""


# "Chunk" of X data, returned by each `Method` above
ChunkX = Union[npt.NDArray[Any], sparse.csr_matrix]


@define
class _SOMAChunk:
    """Return type of ``_ObsAndXSOMAIterator`` that pairs a chunk of ``obs`` rows with the respective rows from the ``X``
    matrix.

    Lifecycle:
        experimental
    """

    obs: pd.DataFrame
    X: ChunkX
    stats: "Stats"

    def __len__(self) -> int:
        return len(self.obs)


Encoders = Dict[str, Encoder]
"""A dictionary of ``Encoder``s keyed by the ``obs`` column name."""


@define
class Stats:
    """Statistics about the data retrieved by ``ExperimentDataPipe`` via SOMA API. This is useful for assessing the read
    throughput of SOMA data.

    Lifecycle:
        experimental
    """

    n_obs: int = 0
    """The total number of obs rows retrieved"""

    nnz: int = 0
    """The total number of values retrieved"""

    elapsed: float = 0
    """The total elapsed time in seconds for retrieving all batches"""

    n_soma_chunks: int = 0
    """The number of chunks retrieved"""

    def __str__(self) -> str:
        return (
            f"{self.n_soma_chunks=}, {self.n_obs=}, {self.nnz=}, "
            f"elapsed={timedelta(seconds=self.elapsed)}"
        )

    def __add__(self, other: "Stats") -> "Stats":
        self.n_obs += other.n_obs
        self.nnz += other.nnz
        self.elapsed += other.elapsed
        self.n_soma_chunks += other.n_soma_chunks
        return self


@define(frozen=True, kw_only=True)
class ExperimentLocator:
    """State required to open the Experiment.

    Necessary as we will likely be invoked across multiple processes.
    """

    uri: str
    tiledb_timestamp_ms: int
    tiledb_config: Dict[str, Union[str, float]]

    @classmethod
    def create(cls, experiment: soma.Experiment) -> "ExperimentLocator":
        return ExperimentLocator(
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


def _tables_to_np(
    tables: Iterator[Tuple[Table, Any]], shape: Tuple[int, int]
) -> typing.Generator[Tuple[npt.NDArray[Any], Any, int], None, None]:
    for tbl, indices in tables:
        row_indices, col_indices, data = (x.to_numpy() for x in tbl.columns)
        nnz = len(data)
        dense_matrix = np.zeros(shape, dtype=data.dtype)
        dense_matrix[row_indices, col_indices] = data
        yield dense_matrix, indices, nnz


class _ObsAndXSOMAIterator(Iterator[_SOMAChunk]):
    """Iterates the SOMA chunks of corresponding ``obs`` and ``X`` data. This is an internal class,
    not intended for public use.
    """

    X: soma.SparseNDArray
    """A handle to the full X data of the SOMA ``Experiment``"""

    obs_joinids_chunks_iter: Iterator[npt.NDArray[np.int64]]

    var_joinids: npt.NDArray[np.int64]
    """The ``var`` joinids to be retrieved from the SOMA ``Experiment``"""

    def __init__(
        self,
        obs: soma.DataFrame,
        X: soma.SparseNDArray,
        obs_column_names: Sequence[str],
        obs_joinids_chunked: List[npt.NDArray[np.int64]],
        var_joinids: npt.NDArray[np.int64],
        shuffle_chunk_count: Optional[int] = None,
        shuffle_rng: Optional[Generator] = None,
        return_sparse_X: bool = False,
    ):
        self.obs = obs
        self.X = X
        self.obs_column_names = obs_column_names
        if shuffle_chunk_count is not None:
            assert shuffle_rng is not None

            # At the start of this step, `obs_joinids_chunked` is a list of one dimensional
            # numpy arrays. Each numpy array corresponds to a chunk of contiguous rows in `obs`.
            # Critically, `obs_joinids_chunked` is randomly ordered where each chunk is
            # from a random section of `obs`.
            # We then take `shuffle_chunk_count` of these in order, concatenate them into
            # a larger numpy array and shuffle this larger numpy array.
            # The result is again a list of numpy arrays.
            self.obs_joinids_chunks_iter = (
                shuffle_rng.permutation(np.concatenate(grouped_chunks))
                for grouped_chunks in list_split(
                    obs_joinids_chunked, shuffle_chunk_count
                )
            )
        else:
            self.obs_joinids_chunks_iter = iter(obs_joinids_chunked)
        self.var_joinids = var_joinids
        self.shuffle_chunk_count = shuffle_chunk_count
        self.return_sparse_X = return_sparse_X

    def __next__(self) -> _SOMAChunk:
        pytorch_logger.debug("Retrieving next SOMA chunk...")
        start_time = time()

        # If no more chunks to iterate through, raise StopIteration, as all iterators do when at end
        obs_joinids_chunk = next(self.obs_joinids_chunks_iter)

        if "soma_joinid" not in self.obs_column_names:
            cols = ["soma_joinid", *self.obs_column_names]
        else:
            cols = list(self.obs_column_names)

        obs_batch = (
            self.obs.read(
                coords=(obs_joinids_chunk,),
                column_names=cols,
            )
            .concat()
            .to_pandas()
            .set_index("soma_joinid")
        )
        assert obs_batch.shape[0] == obs_joinids_chunk.shape[0]

        # handle case of empty result (first batch has 0 rows)
        if len(obs_batch) == 0:
            raise StopIteration

        # reorder obs rows to match obs_joinids_chunk ordering, which may be shuffled
        obs_batch = obs_batch.reindex(obs_joinids_chunk, copy=False)

        # note: the `blockwise` call is employed for its ability to reindex the axes of the sparse matrix,
        # but the blockwise iteration feature is not used (block_size is set to retrieve the chunk as a single block)
        blockwise_iter = self.X.read(
            coords=(obs_joinids_chunk, self.var_joinids)
        ).blockwise(axis=0, size=len(obs_joinids_chunk), eager=False)

        X_batch: ChunkX
        if not self.return_sparse_X:
            res = next(
                _tables_to_np(
                    blockwise_iter.tables(),
                    shape=(obs_batch.shape[0], len(self.var_joinids)),
                )
            )
            X_batch, nnz = res[0], res[2]
        else:
            X_batch = next(blockwise_iter.scipy(compress=True))[0]
            nnz = X_batch.nnz

        assert obs_batch.shape[0] == X_batch.shape[0]

        end_time = time()
        stats = Stats()
        stats.n_obs += X_batch.shape[0]
        stats.nnz += nnz
        stats.elapsed += end_time - start_time
        stats.n_soma_chunks += 1

        pytorch_logger.debug(f"Retrieved SOMA chunk: {stats}")
        return _SOMAChunk(obs=obs_batch, X=X_batch, stats=stats)


def list_split(arr_list: List[Any], sublist_len: int) -> List[List[Any]]:
    """Splits a python list into a list of sublists where each sublist is of size `sublist_len`.
    TODO: Replace with `itertools.batched` when Python 3.12 becomes the minimum supported version.
    """
    i = 0
    result = []
    while i < len(arr_list):
        if (i + sublist_len) >= len(arr_list):
            result.append(arr_list[i:])
        else:
            result.append(arr_list[i : i + sublist_len])

        i += sublist_len

    return result


def run_gc() -> Tuple[Tuple[Any, Any, Any], Tuple[Any, Any, Any], float]:  # noqa: D103
    proc = psutil.Process(os.getpid())

    pre_gc = proc.memory_full_info(), psutil.virtual_memory(), psutil.swap_memory()
    start = time()
    gc.collect()
    gc_elapsed = time() - start
    post_gc = proc.memory_full_info(), psutil.virtual_memory(), psutil.swap_memory()

    pytorch_logger.debug(f"gc:  pre={pre_gc}")
    pytorch_logger.debug(f"gc: post={post_gc}")

    return pre_gc, post_gc, gc_elapsed


class _ObsAndXIterator(Iterator[ObsAndXDatum]):
    """Iterates through a set of ``obs`` and corresponding ``X`` rows, where the rows to be returned are specified by
    the ``obs_tables_iter`` argument. For the specified ``obs` rows, the corresponding ``X`` data is loaded and
    joined together. It is returned from this iterator as 2-tuples of ``X`` and obs Tensors.

    Internally manages the retrieval of data in SOMA-sized chunks, fetching the next chunk of SOMA data as needed.
    Supports fetching the data in an eager manner, where the next SOMA chunk is fetched while the current chunk is
    being read. This is an internal class, not intended for public use.
    """

    soma_chunk_iter: Iterator[_SOMAChunk]
    """The iterator for SOMA chunks of paired obs and X data"""

    soma_chunk: Optional[_SOMAChunk]
    """The current SOMA chunk of obs and X data"""

    i: int = -1
    """Index into current obs ``SOMA`` chunk"""

    def __init__(
        self,
        obs: soma.DataFrame,
        X: soma.SparseNDArray,
        obs_column_names: Sequence[str],
        obs_joinids_chunked: List[npt.NDArray[np.int64]],
        var_joinids: npt.NDArray[np.int64],
        batch_size: int,
        encoders: List[Encoder],
        stats: Stats,
        return_sparse_X: bool,
        use_eager_fetch: bool,
        shuffle_chunk_count: Optional[int] = None,
        shuffle_rng: Optional[Generator] = None,
    ) -> None:
        self.soma_chunk_iter = _ObsAndXSOMAIterator(
            obs,
            X,
            obs_column_names,
            obs_joinids_chunked,
            var_joinids,
            shuffle_chunk_count,
            shuffle_rng,
            return_sparse_X=return_sparse_X,
        )
        if use_eager_fetch:
            self.soma_chunk_iter = _EagerIterator(self.soma_chunk_iter)
        self.soma_chunk = None
        self.var_joinids = var_joinids
        self.batch_size = batch_size
        self.return_sparse_X = return_sparse_X
        self.encoders = encoders
        self.stats = stats
        self.gc_elapsed = 0.0
        self.max_process_mem_usage_bytes = 0
        self.X_dtype = X.schema[2].type.to_pandas_dtype()

    def __next__(self) -> ObsAndXDatum:
        """Read the next torch batch, possibly across multiple soma chunks."""
        obss: list[pd.DataFrame] = []
        Xs: list[ChunkX] = []
        n_obs = 0

        while n_obs < self.batch_size:
            try:
                obs_partial, X_partial = self._read_partial_torch_batch(
                    self.batch_size - n_obs
                )
                n_obs += len(obs_partial)
                obss.append(obs_partial)
                Xs.append(X_partial)
            except StopIteration:
                break

        if len(Xs) == 0:  # If we ran out of data
            raise StopIteration
        else:
            if self.return_sparse_X:
                X = sparse.vstack(Xs)
            else:
                X = np.concatenate(Xs, axis=0)
            obs = pd.concat(obss, axis=0)

        obs_encoded = pd.DataFrame()

        # Add the soma_joinid to the original obs, in case that is requested by the encoders.
        obs["soma_joinid"] = obs.index

        for enc in self.encoders:
            obs_encoded[enc.name] = enc.transform(obs)

        # `to_numpy()` avoids copying the numpy array data
        obs_tensor = torch.from_numpy(obs_encoded.to_numpy())

        if not self.return_sparse_X:
            X_tensor = torch.from_numpy(X)
        else:
            coo = X.tocoo()

            X_tensor = torch.sparse_coo_tensor(
                # Note: The `np.array` seems unnecessary, but PyTorch warns bare array is "extremely slow"
                indices=torch.from_numpy(np.array([coo.row, coo.col])),
                values=coo.data,
                size=coo.shape,
            )

        if self.batch_size == 1:
            X_tensor = X_tensor[0]
            obs_tensor = obs_tensor[0]

        return X_tensor, obs_tensor

    def _read_partial_torch_batch(self, batch_size: int) -> Tuple[pd.DataFrame, ChunkX]:
        """Reads a torch-size batch of data from the current SOMA chunk, returning a torch-size batch whose size may
        contain fewer rows than the requested ``batch_size``. This can happen when the remaining rows in the current
        SOMA chunk are fewer than the requested ``batch_size``.
        """
        if self.soma_chunk is None or not (0 <= self.i < len(self.soma_chunk)):
            # GC memory from previous soma_chunk
            self.soma_chunk = None
            pre_gc, _, gc_elapsed = run_gc()
            self.max_process_mem_usage_bytes = max(
                self.max_process_mem_usage_bytes, pre_gc[0].uss
            )

            self.soma_chunk: _SOMAChunk = next(self.soma_chunk_iter)
            self.stats += self.soma_chunk.stats
            self.gc_elapsed += gc_elapsed
            self.i = 0

            pytorch_logger.debug(
                f"Retrieved SOMA chunk totals: {self.stats}, gc_elapsed={timedelta(seconds=self.gc_elapsed)}"
            )

        obs_batch = self.soma_chunk.obs
        X_chunk = self.soma_chunk.X

        safe_batch_size = min(batch_size, len(obs_batch) - self.i)
        slice_ = slice(self.i, self.i + safe_batch_size)
        assert slice_.stop <= obs_batch.shape[0]

        obs_rows = obs_batch.iloc[slice_]
        assert obs_rows.index.is_unique
        assert safe_batch_size == obs_rows.shape[0]

        X_batch = X_chunk[slice_]

        assert obs_rows.shape[0] == X_batch.shape[0]

        self.i += safe_batch_size

        return obs_rows, X_batch


class ExperimentDataPipe(pipes.IterDataPipe[Dataset[ObsAndXDatum]]):  # type: ignore
    r"""An :class:`torchdata.datapipes.iter.IterDataPipe` that reads ``obs`` and ``X`` data from a
    :class:`tiledbsoma.Experiment`, based upon the specified queries along the ``obs`` and ``var`` axes. Provides an
    iterator over these data when the object is passed to Python's built-in ``iter`` function.

    >>> for batch in iter(ExperimentDataPipe(...)):
            X_batch, y_batch = batch

    The ``batch_size`` parameter controls the number of rows of ``obs`` and ``X`` data that are returned in each
    iteration. If the ``batch_size`` is 1, then each Tensor will have rank 1:

    >>> (tensor([0., 0., 0., 0., 0., 1., 0., 0., 0.]),  # X data
         tensor([2415,    0,    0], dtype=torch.int64)) # obs data, encoded

    For larger ``batch_size`` values, the returned Tensors will have rank 2:

    >>> DataLoader(..., batch_size=3, ...):
        (tensor([[0., 0., 0., 0., 0., 1., 0., 0., 0.],     # X batch
                 [0., 0., 0., 0., 0., 0., 0., 0., 0.],
                 [0., 0., 0., 0., 0., 0., 0., 0., 0.]]),
         tensor([[2415,    0,    0],                       # obs batch
                 [2416,    0,    4],
                 [2417,    0,    3]], dtype=torch.int64))

    The ``return_sparse_X`` parameter controls whether the ``X`` data is returned as a dense or sparse
    :class:`torch.Tensor`. If the model supports use of sparse :class:`torch.Tensor`\ s, this will reduce memory usage.

    The ``obs_column_names`` parameter determines the data columns that are returned in the ``obs`` Tensor. The first
    element is always the ``soma_joinid`` of the ``obs`` :class:`pandas.DataFrame` (or, equivalently, the
    ``soma_dim_0`` of the ``X`` matrix). The remaining elements are the ``obs`` columns specified by
    ``obs_column_names``, and string-typed columns are encoded as integer values. If needed, these values can be decoded
    by obtaining the encoder for a given ``obs`` column name and calling its ``inverse_transform`` method:

    >>> exp_data_pipe.obs_encoders["<obs_attr_name>"].inverse_transform(encoded_values)

    Lifecycle:
        experimental
    """

    _initialized: bool

    _obs_joinids: Optional[npt.NDArray[np.int64]]

    _var_joinids: Optional[npt.NDArray[np.int64]]

    _encoders: List[Encoder]

    _stats: Stats

    _shuffle_rng: Optional[Generator]

    # TODO: Consider adding another convenience method wrapper to construct this object whose signature is more closely
    #  aligned with get_anndata() params (i.e. "exploded" AxisQuery params).
    def __init__(
        self,
        experiment: soma.Experiment,
        measurement_name: str = "RNA",
        X_name: str = "raw",
        obs_query: Optional[soma.AxisQuery] = None,
        var_query: Optional[soma.AxisQuery] = None,
        obs_column_names: Sequence[str] = (),
        batch_size: int = 1,
        shuffle: bool = True,
        seed: Optional[int] = None,
        return_sparse_X: bool = False,
        soma_chunk_size: Optional[int] = 64,
        use_eager_fetch: bool = True,
        encoders: Optional[List[Encoder]] = None,
        shuffle_chunk_count: Optional[int] = 2000,
    ) -> None:
        r"""Construct a new ``ExperimentDataPipe``.

        Args:
            experiment:
                The :class:`tiledbsoma.Experiment` from which to read data.
            measurement_name:
                The name of the :class:`tiledbsoma.Measurement` to read. Defaults to ``"RNA"``.
            X_name:
                The name of the X layer to read. Defaults to ``"raw"``.
            obs_query:
                The query used to filter along the ``obs`` axis. If not specified, all ``obs`` and ``X`` data will
                be returned, which can be very large.
            var_query:
                The query used to filter along the ``var`` axis. If not specified, all ``var`` columns (genes/features)
                will be returned.
            obs_column_names:
                The names of the ``obs`` columns to return. The ``soma_joinid`` index "column" does not need to be
                specified and will always be returned. If not specified, only the ``soma_joinid`` will be returned.
                If custom encoders are passed, this parameter must not be used, since the columns will be inferred
                automatically from the encoders.
            batch_size:
                The number of rows of ``obs`` and ``X`` data to return in each iteration. Defaults to ``1``. A value of
                ``1`` will result in :class:`torch.Tensor` of rank 1 being returns (a single row); larger values will
                result in :class:`torch.Tensor`\ s of rank 2 (multiple rows).
            shuffle:
                Whether to shuffle the ``obs`` and ``X`` data being returned. Defaults to ``True``.
                For performance reasons, shuffling is not performed globally across all rows, but rather in chunks.
                More specifically, we select ``shuffle_chunk_count`` non-contiguous chunks across all the observations
                in the query, concatenate the chunks and shuffle the associated observations.
                The randomness of the shuffling is therefore determined by the
                (``soma_chunk_size``, ``shuffle_chunk_count``) selection. The default values have been determined
                to yield a good trade-off between randomness and performance. Further tuning may be required for
                different type of models. Note that memory usage is correlated to the product
                ``soma_chunk_size * shuffle_chunk_count``.
            seed:
                The random seed used for shuffling. Defaults to ``None`` (no seed). This *must* be specified when using
                :class:`torch.nn.parallel.DistributedDataParallel` to ensure data partitions are disjoint across worker
                processes.
            return_sparse_X:
                Controls whether the ``X`` data is returned as a dense or sparse :class:`torch.Tensor`. As ``X`` data is
                very sparse, setting this to ``True`` will reduce memory usage, if the model supports use of sparse
                :class:`torch.Tensor`\ s. Defaults to ``False``, since sparse :class:`torch.Tensor`\ s are still
                experimental in PyTorch.
            soma_chunk_size:
                The number of ``obs``/``X`` rows to retrieve when reading data from SOMA. This impacts two aspects of
                this class's behavior: 1) The maximum memory utilization, with larger values providing
                better read performance, but also requiring more memory; 2) The granularity of the global shuffling
                step (see ``shuffle`` parameter for details). The default value of 64 works well in conjunction
                with the default ``shuffle_chunk_count`` value.
            use_eager_fetch:
                Fetch the next SOMA chunk of ``obs`` and ``X`` data immediately after a previously fetched SOMA chunk is made
                available for processing via the iterator. This allows network (or filesystem) requests to be made in
                parallel with client-side processing of the SOMA data, potentially improving overall performance at the
                cost of doubling memory utilization. Defaults to ``True``.
            shuffle_chunk_count:
                The number of contiguous blocks (chunks) of rows sampled to then concatenate and shuffle.
                Larger numbers correspond to more randomness per training batch.
                If ``shuffle == False``, this parameter is ignored. Defaults to ``2000``.
            encoders:
                Specify custom encoders to be used. If not specified, a LabelEncoder will be created and
                used for each column in ``obs_column_names``. If specified, only columns for which an encoder
                has been registered will be returned in the ``obs`` tensor. Each encoder needs to have a unique name.
                If this parameter is specified, the ``obs_column_names`` parameter must not be used,
                since the columns will be inferred automatically from the encoders.

        Lifecycle:
            experimental
        """
        self.experiment_locator = ExperimentLocator.create(experiment)
        self.measurement_name = measurement_name
        self.layer_name = X_name
        self.obs_query = obs_query
        self.var_query = var_query
        self.obs_column_names = obs_column_names
        self.batch_size = batch_size
        self.return_sparse_X = return_sparse_X
        self.soma_chunk_size = soma_chunk_size
        self.use_eager_fetch = use_eager_fetch
        self._stats = Stats()
        self._encoders = encoders or []
        self._obs_joinids = None
        self._var_joinids = None
        self._shuffle_chunk_count = shuffle_chunk_count if shuffle else None
        self._shuffle_rng = np.random.default_rng(seed) if shuffle else None
        self._initialized = False
        self.max_process_mem_usage_bytes = 0

        if obs_column_names and encoders:
            raise ValueError(
                "Cannot specify both `obs_column_names` and `encoders`. If `encoders` are specified, columns will be inferred automatically."
            )

        if encoders:
            # Check if names are unique
            if len(encoders) != len({enc.name for enc in encoders}):
                raise ValueError("Encoders must have unique names")

            self.obs_column_names = list(
                dict.fromkeys(itertools.chain(*[enc.columns for enc in encoders]))
            )

    def _init(self) -> None:
        if self._initialized:
            return

        pytorch_logger.debug("Initializing ExperimentDataPipe")

        with self.experiment_locator.open_experiment() as exp:
            query = exp.axis_query(
                measurement_name=self.measurement_name,
                obs_query=self.obs_query,
                var_query=self.var_query,
            )

            # The to_numpy() call is a workaround for a possible bug in TileDB-SOMA:
            # https://github.com/single-cell-data/TileDB-SOMA/issues/1456
            self._obs_joinids = query.obs_joinids().to_numpy()
            self._var_joinids = query.var_joinids().to_numpy()

            self._encoders = self._build_obs_encoders(query)

        self._initialized = True

    @staticmethod
    def _subset_ids_to_partition(
        ids_chunked: List[npt.NDArray[np.int64]],
        partition_index: int,
        num_partitions: int,
    ) -> List[npt.NDArray[np.int64]]:
        """Returns a single partition of the obs_joinids_chunked (a 2D ndarray), based upon the current process's distributed rank and world
        size.
        """
        # subset to a single partition
        # typing does not reflect that is actually a List of 2D NDArrays
        partition_indices = np.array_split(range(len(ids_chunked)), num_partitions)
        partition = [ids_chunked[i] for i in partition_indices[partition_index]]

        if pytorch_logger.isEnabledFor(logging.DEBUG) and len(partition) > 0:
            pytorch_logger.debug(
                f"Process {os.getpid()} handling partition {partition_index + 1} of {num_partitions}, "
                f"partition_size={sum([len(chunk) for chunk in partition])}"
            )

        return partition

    @staticmethod
    def _compute_partitions(
        loader_partition: int,
        loader_partitions: int,
        dist_partition: int,
        num_dist_partitions: int,
    ) -> Tuple[int, int]:
        # NOTE: Can alternately use a `worker_init_fn` to split among workers split workload
        total_partitions = num_dist_partitions * loader_partitions
        partition = dist_partition * loader_partitions + loader_partition
        return partition, total_partitions

    def __iter__(self) -> Iterator[ObsAndXDatum]:
        self._init()
        assert self._obs_joinids is not None
        assert self._var_joinids is not None

        if self.soma_chunk_size is None:
            # set soma_chunk_size to utilize ~1 GiB of RAM per SOMA chunk; assumes 95% X data sparsity, 8 bytes for the
            # X value and 8 bytes for the sparse matrix indices, and a 100% working memory overhead (2x).
            X_row_memory_size = 0.05 * len(self._var_joinids) * 8 * 3 * 2
            self.soma_chunk_size = int((1 * 1024**3) / X_row_memory_size)
        pytorch_logger.debug(f"Using {self.soma_chunk_size=}")

        if (
            self.return_sparse_X
            and torch.utils.data.get_worker_info()
            and torch.utils.data.get_worker_info().num_workers > 0
        ):
            raise NotImplementedError(
                "torch does not work with sparse tensors in multi-processing mode "
                "(see https://github.com/pytorch/pytorch/issues/20248)"
            )

        # chunk the obs joinids into batches of size soma_chunk_size
        obs_joinids_chunked = self._chunk_ids(self._obs_joinids, self.soma_chunk_size)

        # globally shuffle the chunks, if requested
        if self._shuffle_rng:
            self._shuffle_rng.shuffle(obs_joinids_chunked)

        # subset to a single partition, as needed for distributed training and multi-processing datat loading
        worker_info = torch.utils.data.get_worker_info()
        partition, partitions = self._compute_partitions(
            loader_partition=worker_info.id if worker_info else 0,
            loader_partitions=worker_info.num_workers if worker_info else 1,
            dist_partition=dist.get_rank() if dist.is_initialized() else 0,
            num_dist_partitions=dist.get_world_size() if dist.is_initialized() else 1,
        )
        obs_joinids_chunked_partition: List[
            npt.NDArray[np.int64]
        ] = self._subset_ids_to_partition(obs_joinids_chunked, partition, partitions)

        with self.experiment_locator.open_experiment() as exp:
            X = exp.ms[self.measurement_name].X[self.layer_name]
            if not isinstance(X, soma.SparseNDArray):
                raise NotImplementedError(
                    "ExperimentDataPipe only supported on X layers which are of type SparseNDArray"
                )

            obs_and_x_iter = _ObsAndXIterator(
                obs=exp.obs,
                X=X,
                obs_column_names=self.obs_column_names,
                obs_joinids_chunked=obs_joinids_chunked_partition,
                var_joinids=self._var_joinids,
                batch_size=self.batch_size,
                encoders=self._encoders,
                stats=self._stats,
                return_sparse_X=self.return_sparse_X,
                use_eager_fetch=self.use_eager_fetch,
                shuffle_rng=self._shuffle_rng,
                shuffle_chunk_count=self._shuffle_chunk_count,
            )

            yield from obs_and_x_iter

            self.max_process_mem_usage_bytes = (
                obs_and_x_iter.max_process_mem_usage_bytes
            )
            pytorch_logger.debug(
                "max process memory usage="
                f"{self.max_process_mem_usage_bytes / (1024 ** 3):.3f} GiB"
            )

    @staticmethod
    def _chunk_ids(
        ids: npt.NDArray[np.int64], chunk_size: int
    ) -> List[npt.NDArray[np.int64]]:
        num_chunks = max(1, ceil(len(ids) / chunk_size))
        pytorch_logger.debug(
            f"Shuffling {len(ids)} obs joinids into {num_chunks} chunks of {chunk_size}"
        )
        return np.array_split(ids, num_chunks)

    def __len__(self) -> int:
        self._init()
        assert self._obs_joinids is not None

        div, rem = divmod(len(self._obs_joinids), self.batch_size)
        return div + bool(rem)

    def __getitem__(self, index: int) -> ObsAndXDatum:
        raise NotImplementedError("IterDataPipe can only be iterated")

    def _build_obs_encoders(self, query: soma.ExperimentAxisQuery) -> List[Encoder]:
        pytorch_logger.debug("Initializing encoders")

        encoders = []

        if "soma_joinid" not in self.obs_column_names:
            cols = ["soma_joinid", *self.obs_column_names]
        else:
            cols = list(self.obs_column_names)

        obs = query.obs(column_names=cols).concat().to_pandas()

        if self._encoders:
            # Fit all the custom encoders with obs
            for enc in self._encoders:
                enc.fit(obs)
                encoders.append(enc)
        else:
            # Create one LabelEncoder for each column, and fit it with obs
            for col in self.obs_column_names:
                enc = LabelEncoder(col)
                enc.fit(obs)
                encoders.append(enc)

        return encoders

    # TODO: This does not work in multiprocessing mode, as child process's stats are not collected
    def stats(self) -> Stats:
        """Get data loading stats for this :class:`tiledbsoma.ml.pytorch.ExperimentDataPipe`.

        Returns:
            The :class:`tiledbsoma.ml.pytorch.Stats` object for this
            :class:`tiledbsoma.ml.pytorch.ExperimentDataPipe`.

        Lifecycle:
            experimental
        """
        return self._stats

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
        self._init()
        assert self._obs_joinids is not None
        assert self._var_joinids is not None

        return len(self._obs_joinids), len(self._var_joinids)

    @property
    def obs_encoders(self) -> Encoders:
        """Returns a dictionary of :class:`sklearn.preprocessing.LabelEncoder` objects, keyed on ``obs`` column names,
        which were used to encode the ``obs`` column values.

        These encoders can be used to decode the encoded values as follows:

        >>> exp_data_pipe.obs_encoders["<obs_attr_name>"].inverse_transform(encoded_values)

        Returns:
            A ``Dict[str, LabelEncoder]``, mapping column names to :class:`sklearn.preprocessing.LabelEncoder` objects.
        """
        self._init()
        assert self._encoders is not None

        return {enc.name: enc for enc in self._encoders}


# Note: must be a top-level function (and not a lambda), to play nice with multiprocessing pickling
def _collate_noop(x: Any) -> Any:
    return x


# TODO: Move into somacore.ExperimentAxisQuery
def experiment_dataloader(
    datapipe: pipes.IterDataPipe,
    num_workers: int = 0,
    **dataloader_kwargs: Any,
) -> DataLoader:
    """Factory method for :class:`torch.utils.data.DataLoader`. This method can be used to safely instantiate a
    :class:`torch.utils.data.DataLoader` that works with :class:`tiledbsoma.ml.pytorch.ExperimentDataPipe`,
    since some of the :class:`torch.utils.data.DataLoader` constructor parameters are not applicable when using a
    :class:`torchdata.datapipes.iter.IterDataPipe` (``shuffle``, ``batch_size``, ``sampler``, ``batch_sampler``,
    ``collate_fn``).

    Args:
        datapipe:
            An :class:`torchdata.datapipes.iter.IterDataPipe`, which can be an
            :class:`tiledbsoma.ml.pytorch.ExperimentDataPipe` or any other
            :class:`torchdata.datapipes.iter.IterDataPipe` that has been chained to the
            :class:`tiledbsoma.ml.pytorch.ExperimentDataPipe`.
        num_workers:
            Number of worker processes to use for data loading. If ``0``, data will be loaded in the main process.
        **dataloader_kwargs:
            Additional keyword arguments to pass to the :class:`torch.utils.data.DataLoader` constructor,
            except for ``shuffle``, ``batch_size``, ``sampler``, ``batch_sampler``, and ``collate_fn``, which are not
            supported when using :class:`tiledbsoma.ml.pytorch.ExperimentDataPipe`.

    Returns:
        A :class:`torch.utils.data.DataLoader`.

    Raises:
        ValueError: if any of the ``shuffle``, ``batch_size``, ``sampler``, ``batch_sampler``, or ``collate_fn`` params
            are passed as keyword arguments.

    Lifecycle:
        experimental
    """
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

    return DataLoader(
        datapipe,
        batch_size=None,  # batching is handled by our ExperimentDataPipe
        num_workers=num_workers,
        # avoid use of default collator, which adds an extra (3rd) dimension to the tensor batches
        collate_fn=_collate_noop,
        # shuffling is handled by our ExperimentDataPipe
        shuffle=False,
        **dataloader_kwargs,
    )


def _init_multiprocessing() -> None:
    """Ensures use of "spawn" for starting child processes with multiprocessing.

    Forked processes are known to be problematic:
      https://pytorch.org/docs/stable/notes/multiprocessing.html#avoiding-and-fighting-deadlocks
    Also, CUDA does not support forked child processes:
      https://pytorch.org/docs/stable/notes/multiprocessing.html#cuda-in-multiprocessing

    """
    torch.multiprocessing.set_start_method("fork", force=True)
    orig_start_method = torch.multiprocessing.get_start_method()
    if orig_start_method != "spawn":
        if orig_start_method:
            pytorch_logger.warning(
                "switching torch multiprocessing start method from "
                f'"{torch.multiprocessing.get_start_method()}" to "spawn"'
            )
        torch.multiprocessing.set_start_method("spawn", force=True)
