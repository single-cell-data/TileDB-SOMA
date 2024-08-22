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
    List,
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
import psutil
import pyarrow as pa
import scipy.sparse as sparse
import torch
import torchdata
from somacore.query._eager_iter import EagerIterator as _EagerIterator

import tiledbsoma as soma

from .encoders import Encoder, LabelEncoder

pytorch_logger = logging.getLogger("tiledbsoma.ml.pytorch")

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


Encoders = Dict[str, Encoder]
"""A dictionary of ``Encoder``s keyed by the ``obs`` column name."""


if TYPE_CHECKING:
    PsUtilMemInfo = Tuple[psutil.pfullmem, psutil.svmem, psutil.sswap]
else:
    PsUtilMemInfo = Tuple[Any]


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
        obs_column_names: Sequence[str] = (),
        batch_size: int = 1,
        shuffle: bool = True,
        io_batch_size: int = 2**17,
        shuffle_chunk_size: int = 64,
        seed: int | None = None,
        use_eager_fetch: bool = True,
        encoders: List[Encoder] | None = None,
        partition: bool = True,
    ):
        super().__init__()

        # Anything set in the instance needs to be picklable for multi-process users
        self.experiment_locator = _ExperimentLocator.create(experiment)
        self.measurement_name = measurement_name
        self.layer_name = X_name
        self.obs_query = obs_query
        self.var_query = var_query
        self.obs_column_names = obs_column_names
        self.batch_size = batch_size
        self.io_batch_size = io_batch_size
        self.shuffle = shuffle
        self.use_eager_fetch = use_eager_fetch
        # XXX - TODO: when/if we add X encoders, how will they be differentiated from obs encoders?  Naming?
        self._encoders = encoders or []
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

        if pytorch_logger.isEnabledFor(logging.DEBUG) and self.partition:
            pytorch_logger.debug(
                f"Process {os.getpid()} handling partition {partition + 1} of {num_partitions}, "
                f"partition_size={sum([len(chunk) for chunk in obs_partition_joinids])}"
            )

        return obs_joinid_iter

    def _init_once(self) -> None:
        """One-time per worker initialization. All operations be idempotent in order to support pipe reset()."""
        if self._initialized:
            return

        pytorch_logger.debug(
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
            self._encoders = self._build_encoders(query)

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
            _mini_batch_iter = self._encoded_mini_batch_iter(
                exp.obs, X, obs_joinid_iter
            )
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
            # obs_coords = np.concatenate(obs_coord_chunks)
            obs_shuffled_coords = (
                obs_coords
                if self._shuffle_rng is None
                else self._shuffle_rng.permuted(obs_coords)
            )
            pytorch_logger.debug(
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

            X_tbl_iter = X.read(coords=(obs_coords, self._var_joinids)).tables()
            X_tbl = pa.concat_tables(X_tbl_iter)
            obs_indexer = soma.IntIndexer(obs_shuffled_coords, context=X.context)
            X_tbl = pa.Table.from_pydict(
                {
                    "soma_dim_0": obs_indexer.get_indexer(
                        X_tbl["soma_dim_0"].to_numpy()
                    ),
                    "soma_dim_1": var_indexer.get_indexer(
                        X_tbl["soma_dim_1"].to_numpy()
                    ),
                    "soma_data": X_tbl["soma_data"].to_numpy(),
                }
            )
            X_io_batch = sparse.csr_array(
                _D_IJ(X_tbl),
                shape=(len(obs_coords), len(self._var_joinids)),
                copy=False,
            )

            del X_tbl, X_tbl_iter
            del obs_coords, obs_shuffled_coords, obs_indexer
            _run_gc()
            tm = time.perf_counter() - st_time
            pytorch_logger.debug(
                f"Retrieved SOMA IO batch, took {tm:.2f}sec, {X_io_batch.shape[0]/tm:0.1f} samples/sec"
            )
            yield X_io_batch, obs_io_batch

    def _mini_batch_iter(
        self,
        obs: soma.DataFrame,
        X: soma.SparseNDArray,
        obs_joinid_iter: Iterator[npt.NDArray[np.int64]],
    ) -> Iterator[Tuple[sparse.csr_array | sparse.csr_matrix, pd.DataFrame]]:
        """Break IO batches into shuffled mini-batch-sized chunks, still in internal format (CSR, Pandas).

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
                    yield result
                    result = None

        else:
            # yield a remnant, if any
            if result is not None:
                yield result

    def _encoded_mini_batch_iter(
        self,
        obs: soma.DataFrame,
        X: soma.SparseNDArray,
        obs_joinid_iter: Iterator[npt.NDArray[np.int64]],
    ) -> Iterator[XObsDatum]:
        """Apply encoding on top of the mini batches.

        Returns numpy encodings (X, obs).
        """
        _encoded_mini_batch_iter = self._mini_batch_iter(obs, X, obs_joinid_iter)
        if self.use_eager_fetch:
            _encoded_mini_batch_iter = _EagerIterator(
                _encoded_mini_batch_iter, pool=X.context.threadpool
            )

        for X_mini_batch, obs_mini_batch in _encoded_mini_batch_iter:
            # TODO - X encoding
            X_encoded = _csr_to_dense(
                X_mini_batch
            )  # same as X_mini_batch.todense(), which is SLOW

            # Obs encoding
            obs_encoded = (
                pd.DataFrame(
                    {enc.name: enc.transform(obs_mini_batch) for enc in self._encoders}
                )
                if self._encoders
                else obs_mini_batch
            )

            del obs_mini_batch, X_mini_batch
            yield X_encoded, obs_encoded

    def _build_encoders(
        self, query: soma.ExperimentAxisQuery[soma.Experiment]  # type: ignore[type-var]
    ) -> List[Encoder]:
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

    @property
    def encoders(self) -> Encoders:
        """Returns a dictionary of :class:`sklearn.preprocessing.LabelEncoder` objects, keyed on ``obs`` column names,
        which were used to encode the ``obs`` column values.

        These encoders can be used to decode the encoded values as follows:

        >>> exp_data_pipe.encoders["<obs_attr_name>"].inverse_transform(encoded_values)

        Returns:
            A ``Dict[str, LabelEncoder]``, mapping column names to :class:`sklearn.preprocessing.LabelEncoder` objects.
        """
        self._init_once()
        assert self._encoders is not None
        return {enc.name: enc for enc in self._encoders}


class ExperimentAxisQueryDataPipe(
    torchdata.datapipes.iter.IterDataPipe[  # type:ignore[misc]
        torch.utils.data.dataset.Dataset[XObsTensorDatum]
    ],
):
    """An :class:`torchdata.datapipes.iter.IterDataPipe` which reads ``obs`` and ``X`` data from a
    :class:`tiledbsoma.Experiment`, based upon the specified queries along the ``obs`` and ``var`` axes. Provides
    a standard Python iterable interface.

    >>> for batch in ExperimentAxisQueryDataPipe(...):
            X_batch, obs_batch = batch

    **WARNING:** :class:`torchdata.datapipes` is deprecated as of version 0.9 (July 2024), and is slated for removal
    in a future release (late 2024). It is recommended that new code utilize :class:`ExperimentAxisQueryIterableDataset`.
    Older code should pin the torchdata version to 0.9 or older. For more information, see
    https://github.com/pytorch/data/issues/1196.

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

    The ``obs_column_names`` parameter determines the data columns that are returned in the ``obs`` Tensor.
    User-specified encoders may be provided - when not provided, the ``X`` batch will not be encoded, and the
    ``obs`` batch will be encoded with a simple label encoder.

    If needed, encoded valeus can be decoded by calling ``inverse_transform`` method on the encoder, e.g.,

    >>> exp_data_pipe.encoders["<obs_attr_name>"].inverse_transform(encoded_values)

    Lifecycle:
        experimental
    """

    def __init__(
        self,
        experiment: soma.Experiment,
        measurement_name: str = "RNA",
        X_name: str = "raw",
        obs_query: soma.AxisQuery | None = None,
        var_query: soma.AxisQuery | None = None,
        obs_column_names: Sequence[str] = (),
        batch_size: int = 1,
        shuffle: bool = True,
        seed: int | None = None,
        io_batch_size: int = 2**17,
        shuffle_chunk_size: int = 64,
        use_eager_fetch: bool = True,
        encoders: List[Encoder] | None = None,
    ):
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
            encoders=encoders,
            shuffle_chunk_size=shuffle_chunk_size,
            # ---
            partition=True,
        )

    def __iter__(self) -> Iterator[XObsTensorDatum]:
        batch_size = self._exp_iter.batch_size
        for X, obs in self._exp_iter:
            X_tensor, obs_tensor = torch.from_numpy(X), torch.from_numpy(obs.to_numpy())
            if batch_size == 1:
                X_tensor = X_tensor[0]
                obs_tensor = obs_tensor[0]

            yield X_tensor, obs_tensor

    def __len__(self) -> int:
        return self._exp_iter.__len__()

    @property
    def shape(self) -> Tuple[int, int]:
        return self._exp_iter.shape

    @property
    def encoders(self) -> Encoders:
        return self._exp_iter.encoders


class ExperimentAxisQueryIterableDataset(
    torch.utils.data.IterableDataset[XObsNpDatum]  # type:ignore[misc]
):
    def __init__(
        self,
        experiment: soma.Experiment,
        measurement_name: str = "RNA",
        X_name: str = "raw",
        obs_query: soma.AxisQuery | None = None,
        var_query: soma.AxisQuery | None = None,
        obs_column_names: Sequence[str] = (),
        batch_size: int = 1,  # XXX add docstring noting values >1 will not work with default collator
        shuffle: bool = True,
        seed: int | None = None,
        io_batch_size: int = 2**17,
        shuffle_chunk_size: int = 64,
        use_eager_fetch: bool = True,
        encoders: List[Encoder] | None = None,
    ):
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
            encoders=encoders,
            shuffle_chunk_size=shuffle_chunk_size,
            # ---
            partition=True,
        )

    def __iter__(self) -> Iterator[XObsNpDatum]:
        batch_size = self._exp_iter.batch_size
        for X, obs in self._exp_iter:
            obs_np: NDArrayNumber = obs.to_numpy()
            if batch_size == 1:
                X = X[0]
                obs_np = obs_np[0]

            yield X, obs_np

    def __len__(self) -> int:
        return self._exp_iter.__len__()

    @property
    def shape(self) -> Tuple[int, int]:
        return self._exp_iter.shape

    @property
    def encoders(self) -> Encoders:
        return self._exp_iter.encoders


def _collate_ndarray_to_tensor(
    datum: Sequence[NDArrayNumber | torch.Tensor],
) -> Tuple[torch.Tensor, ...]:
    """Default torch.utils.data.DataLoader collate_fn for ``experiment_dataloader`` -- converts ndarray to a Tensor.

    Must be a top-level function to play nice with picking for multiprocessing use cases.
    """
    return tuple(
        torch.from_numpy(d) if not isinstance(d, torch.Tensor) else d for d in datum
    )


def experiment_dataloader(
    ds: torchdata.datapipes.iter.IterDataPipe | torch.utils.data.IterableDataset,
    num_workers: int = 0,
    **dataloader_kwargs: Any,
) -> torch.utils.data.DataLoader:
    """Factory method for :class:`torch.utils.data.DataLoader`. This method can be used to safely instantiate a
    :class:`torch.utils.data.DataLoader` that works with :class:`tiledbsoma.ml.ExperimentAxisQueryIterableDataset`.

    Several :class:`torch.utils.data.DataLoader` parameters are disallowed as they are either non-performant and better
    specified on the :class:`tiledbsoma.ml.ExperimentAxisQueryIterableDataset` (``shuffle``, ``batch_size``) or are not compatible
    (``sampler``, ``batch_sampler``, ``collate_fn``).

    Args:
        datapipe:
            An :class:`torch.util.data.IterableDataset`, which can be an
            :class:`tiledbsoma.ml.ExperimentAxisQueryIterableDataset`,
            :class:`tiledbsoma.ml.ExperimentAxisQueryIterDataPipe` or any other
            :class:`torch.util.data.IterableDataset` that has been chained to the
            :class:`tiledbsoma.ml.ExperimentAxisQueryIterableDataset`.
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

    return torch.utils.data.DataLoader(
        ds,
        batch_size=None,  # batching is handled by our ExperimentAxisQueryIterableDataset
        num_workers=num_workers,
        collate_fn=_collate_ndarray_to_tensor,
        # shuffling is handled by ExperimentAxisQueryIterableDataset
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


def _D_IJ(
    tbl: pa.Table,
) -> Tuple[NDArrayNumber, Tuple[npt.NDArray[np.int64], npt.NDArray[np.int64]]]:
    """Given SOMA-style Pyarrow Table of COO sparse array data, return tuple (D, (I, J)) vectors."""
    d = tbl["soma_data"].to_numpy()
    i = tbl["soma_dim_0"].to_numpy()
    j = tbl["soma_dim_1"].to_numpy()
    return d, (i, j)


if sys.version_info >= (3, 12):
    _batched = itertools.batched

else:

    def _batched(
        iterable: Iterable[_T_co], n: int, *, strict: bool = False
    ) -> Iterator[Tuple[_T_co, ...]]:
        """Same as the Python 3.12+ itertools.batched, but polyfilled for old implementations."""
        if n < 1:
            raise ValueError("n must be at least one")
        it = iter(iterable)
        while batch := tuple(islice(it, n)):
            if strict and len(batch) != n:
                raise ValueError("batched(): incomplete batch")
            yield batch


def _run_gc() -> Tuple[PsUtilMemInfo, PsUtilMemInfo, float]:
    """Run Python GC and log stats."""
    proc = psutil.Process(os.getpid())

    pre_gc = proc.memory_full_info(), psutil.virtual_memory(), psutil.swap_memory()
    start = time.time()
    gc.collect()
    gc_elapsed = time.time() - start
    post_gc = proc.memory_full_info(), psutil.virtual_memory(), psutil.swap_memory()

    pytorch_logger.debug(f"gc:  pre={pre_gc}")
    pytorch_logger.debug(f"gc: post={post_gc}")

    return pre_gc, post_gc, gc_elapsed


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
            pytorch_logger.warning(
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
    assert isinstance(sp, (sparse.csr_array, sparse.csr_matrix))
    return cast(
        NDArrayNumber,
        _csr_to_dense_inner(
            sp.indptr, sp.indices, sp.data, np.zeros(sp.shape, dtype=sp.dtype)
        ),
    )
