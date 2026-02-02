# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.
from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Mapping, Sequence
from typing import Any, Union

import numpy as np
import numpy.typing as npt
import pyarrow as pa
from anndata import AnnData
from typing_extensions import Protocol, Self, TypedDict

from tiledbsoma._core import DataFrame, ReadIter, SparseRead, measurement
from tiledbsoma._core import types as base_types
from tiledbsoma._core.options import BatchSize, PlatformConfig, ReadPartitions, ResultOrder, ResultOrderStr

_RO_AUTO = ResultOrder.AUTO
_UNBATCHED = BatchSize()


class AxisColumnNames(TypedDict, total=False):
    """Specifies column names for experiment axis query read operations.

    Lifecycle: maturing
    """

    obs: Sequence[str] | None
    """obs columns to use. All columns if ``None`` or not present."""
    var: Sequence[str] | None
    """var columns to use. All columns if ``None`` or not present."""


class ExperimentAxisQuery(ABC):
    """Axis-based query against a SOMA Experiment.

    ExperimentAxisQuery allows easy selection and extraction of data from a
    single :class:`Measurement` in an :class:`Experiment`, by obs/var (axis) coordinates
    and/or value filter.

    The primary use for this class is slicing :class:`Experiment` ``X`` layers by obs or
    var value and/or coordinates. Slicing on :class:`SparseNDArray` ``X`` matrices is
    supported; :class:`DenseNDArray` is not supported at this time.

    Lifecycle: maturing
    """

    @abstractmethod
    def obs(
        self,
        *,
        column_names: Sequence[str] | None = None,
        batch_size: BatchSize = _UNBATCHED,
        partitions: ReadPartitions | None = None,
        result_order: ResultOrderStr = _RO_AUTO,
        platform_config: PlatformConfig | None = None,
    ) -> ReadIter[pa.Table]:
        """Returns ``obs`` as an `Arrow table
        <https://arrow.apache.org/docs/python/generated/pyarrow.Table.html>`_
        iterator.

        Lifecycle: maturing
        """
        ...

    @abstractmethod
    def var(
        self,
        *,
        column_names: Sequence[str] | None = None,
        batch_size: BatchSize = _UNBATCHED,
        partitions: ReadPartitions | None = None,
        result_order: ResultOrderStr = _RO_AUTO,
        platform_config: PlatformConfig | None = None,
    ) -> ReadIter[pa.Table]:
        """Returns ``var`` as an `Arrow table
        <https://arrow.apache.org/docs/python/generated/pyarrow.Table.html>`_
        iterator.

        Lifecycle: maturing
        """
        ...

    @abstractmethod
    def obs_joinids(self) -> pa.IntegerArray:
        """Returns ``obs`` ``soma_joinids`` as an Arrow array.

        Lifecycle: maturing
        """
        ...

    @abstractmethod
    def var_joinids(self) -> pa.IntegerArray:
        """Returns ``var`` ``soma_joinids`` as an Arrow array.

        Lifecycle: maturing
        """
        ...

    @property
    @abstractmethod
    def n_obs(self) -> int:
        """The number of ``obs`` axis query results.

        Lifecycle: maturing
        """
        ...

    @property
    @abstractmethod
    def n_vars(self) -> int:
        """The number of ``var`` axis query results.

        Lifecycle: maturing
        """
        ...

    @property
    @abstractmethod
    def indexer(self) -> AxisIndexer:
        """A ``soma_joinid`` indexer for both ``obs`` and ``var`` axes.

        Lifecycle: maturing
        """
        ...

    @abstractmethod
    def X(
        self,
        layer_name: str,
        *,
        batch_size: BatchSize = _UNBATCHED,
        partitions: ReadPartitions | None = None,
        result_order: ResultOrderStr = _RO_AUTO,
        platform_config: PlatformConfig | None = None,
    ) -> SparseRead:
        """Returns an ``X`` layer as a sparse read.

        Args:
            layer_name: The X layer name to return.
            batch_size: The size of batches that should be returned from a read.
                See :class:`BatchSize` for details.
            partitions: Specifies that this is part of a partitioned read,
                and which partition to include, if present.
            result_order: the order to return results, specified as a
                :class:`~ResultOrder` or its string value.
            platform_config: platform-specific configuration; keys are SOMA
                implementation names.

        Lifecycle: maturing
        """
        ...

    @abstractmethod
    def obsp(self, layer: str) -> SparseRead:
        """Returns an ``obsp`` layer as a sparse read.

        Lifecycle: maturing
        """
        ...

    @abstractmethod
    def varp(self, layer: str) -> SparseRead:
        """Returns a ``varp`` layer as a sparse read.

        Lifecycle: maturing
        """
        ...

    @abstractmethod
    def obsm(self, layer: str) -> SparseRead:
        """Returns an ``obsm`` layer as a sparse read.

        Lifecycle: maturing
        """
        ...

    @abstractmethod
    def varm(self, layer: str) -> SparseRead:
        """Returns a ``varm`` layer as a sparse read.

        Lifecycle: maturing
        """
        ...

    @abstractmethod
    def obs_scene_ids(self) -> pa.Array:
        """Returns a pyarrow array with scene ids that contain obs from this
        query.

        Lifecycle: experimental
        """
        ...

    @abstractmethod
    def var_scene_ids(self) -> pa.Array:
        """Return a pyarrow array with scene ids that contain var from this
        query.

        Lifecycle: experimental
        """
        ...

    @abstractmethod
    def to_anndata(
        self,
        X_name: str,
        *,
        column_names: AxisColumnNames | None = None,
        X_layers: Sequence[str] = (),
        obsm_layers: Sequence[str] = (),
        obsp_layers: Sequence[str] = (),
        varm_layers: Sequence[str] = (),
        varp_layers: Sequence[str] = (),
        drop_levels: bool = False,
    ) -> AnnData:
        """Executes the query and return result as an ``AnnData`` in-memory object.

        Args:
            X_name: The X layer to read and return in the ``X`` slot.
            column_names: The columns in the ``var`` and ``obs`` dataframes
                to read.
            X_layers: Additional X layers to read and return
                in the ``layers`` slot.
            obsm_layers:
                Additional obsm layers to read and return in the obsm slot.
            obsp_layers:
                Additional obsp layers to read and return in the obsp slot.
            varm_layers:
                Additional varm layers to read and return in the varm slot.
            varp_layers:
                Additional varp layers to read and return in the varp slot.
            drop_levels:
                Indicate whether unused categories on axis frames should be
                dropped. By default, False, the categories which are present
                in the SOMA Experiment and not present in the query output
                are not dropped.

        Lifecycle: maturing
        """
        ...

    # Context management

    @abstractmethod
    def close(self) -> None:
        """Releases resources associated with this query.

        This method must be idempotent.

        Lifecycle: maturing
        """
        ...

    @abstractmethod
    def __enter__(self) -> Self: ...

    @abstractmethod
    def __exit__(self, *_: Any) -> None: ...  # noqa: ANN401


Numpyable = Union[pa.Array, pa.ChunkedArray, npt.NDArray[np.int64]]
"""Things that can be converted to a NumPy array."""


class AxisIndexer(ABC):
    """Given a query, provides index-building services for obs/var axis.

    Lifecycle: maturing
    """

    @abstractmethod
    def by_obs(self, coords: Numpyable) -> npt.NDArray[np.intp]:
        """Reindex the coords (soma_joinids) over the ``obs`` axis."""
        ...

    @abstractmethod
    def by_var(self, coords: Numpyable) -> npt.NDArray[np.intp]:
        """Reindex for the coords (soma_joinids) over the ``var`` axis."""
        ...


class Experimentish(Protocol):
    """The API we need from an Experiment."""

    @property
    def ms(self) -> Mapping[str, measurement.Measurement]: ...  # type: ignore[type-arg]

    @property
    def obs(self) -> DataFrame: ...

    @property
    def context(self) -> base_types.ContextBase | None: ...

    @property
    def obs_spatial_presence(self) -> DataFrame: ...
