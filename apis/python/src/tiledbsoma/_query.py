# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""Implementation of a SOMA Experiment."""

from __future__ import annotations

import enum
import json
import warnings
from collections.abc import Sequence
from concurrent.futures import Future, ThreadPoolExecutor
from threading import Lock
from typing import TYPE_CHECKING, Any, Callable, Literal, Protocol, TypeVar, Union

import attrs
import numpy as np
import numpy.typing as npt
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pacomp
import scipy.sparse as sp
from anndata import AnnData
from typing_extensions import Self, TypedDict

from ._axis import AxisQuery
from ._collection import Collection
from ._constants import SOMA_DATAFRAME_ORIGINAL_INDEX_NAME_JSON, SPATIAL_DISCLAIMER
from ._core_iters import _RO_AUTO, ReadIter, SparseRead
from ._core_options import (
    BatchSize,
    PlatformConfig,
    ReadPartitions,
    ResultOrder,
    ResultOrderStr,
)
from ._dask.load import SOMADaskConfig, load_daskarray
from ._dataframe import DataFrame
from ._exception import SOMAError

if TYPE_CHECKING:
    from ._experiment import Experiment

from ._fastercsx import CompressedMatrix
from ._measurement import Measurement
from ._sparse_nd_array import SparseNDArray
from ._types import IndexFactory, IndexLike
from ._util import MISSING, Sentinel, _df_set_index, _resolve_futures

_T = TypeVar("_T")
_T_co = TypeVar("_T_co", covariant=True)

Numpyable = Union[pa.Array, pa.ChunkedArray, npt.NDArray[np.int64]]
"""Things that can be converted to a NumPy array."""


class AxisColumnNames(TypedDict, total=False):
    """Specifies column names for experiment axis query read operations.

    Lifecycle: maturing
    """

    obs: Sequence[str] | None
    """obs columns to use. All columns if ``None`` or not present."""
    var: Sequence[str] | None
    """var columns to use. All columns if ``None`` or not present."""


class _HasObsVar(Protocol[_T_co]):
    """Something which has an ``obs`` and ``var`` field.

    Used to give nicer type inference in :meth:`Axis.getattr_from`.
    """

    @property
    def obs(self) -> _T_co: ...

    @property
    def var(self) -> _T_co: ...


class AxisName(enum.Enum):
    OBS = "obs"
    VAR = "var"

    @property
    def value(self) -> Literal["obs", "var"]:
        return super().value


@attrs.define
class AxisIndexer:
    """Given a query, provides index-building services for obs/var axis.

    Lifecycle: maturing
    """

    query: ExperimentAxisQuery
    _index_factory: IndexFactory
    _cached_obs: IndexLike | None = None
    _cached_var: IndexLike | None = None
    _obs_lock: Lock = attrs.field(factory=Lock)
    _var_lock: Lock = attrs.field(factory=Lock)

    @property
    def _obs_index(self) -> IndexLike:
        """Private. Return an index for the ``obs`` axis."""
        with self._obs_lock:
            if self._cached_obs is None:
                self._cached_obs = self._index_factory(self.query.obs_joinids().to_numpy())
        return self._cached_obs

    @property
    def _var_index(self) -> IndexLike:
        """Private. Return an index for the ``var`` axis."""
        with self._var_lock:
            if self._cached_var is None:
                self._cached_var = self._index_factory(self.query.var_joinids().to_numpy())
        return self._cached_var

    def by_obs(self, coords: Numpyable) -> npt.NDArray[np.intp]:
        """Reindex the coords (soma_joinids) over the ``obs`` axis."""
        return self._obs_index.get_indexer(_to_numpy(coords))

    def by_var(self, coords: Numpyable) -> npt.NDArray[np.intp]:
        """Reindex for the coords (soma_joinids) over the ``var`` axis."""
        return self._var_index.get_indexer(_to_numpy(coords))


def _to_numpy(it: Numpyable) -> npt.NDArray[np.int64]:
    if isinstance(it, np.ndarray):
        return it
    return it.to_numpy()


class ExperimentAxisQuery:
    """Axis-based query against a SOMA Experiment.

    ExperimentAxisQuery allows easy selection and extraction of data from a
    single :class:`Measurement` in an :class:`Experiment`, by obs/var (axis) coordinates
    and/or value filter.

    The primary use for this class is slicing :class:`Experiment` ``X`` layers by obs or
    var value and/or coordinates. Slicing on :class:`SparseNDArray` ``X`` matrices is
    supported; :class:`DenseNDArray` is not supported at this time.

    IMPORTANT: this class is not thread-safe.

    IMPORTANT: this query class assumes it can store the full result of both
    axis dataframe queries in memory, and only provides incremental access to
    the underlying X NDArray. API features such as ``n_obs`` and ``n_vars``
    codify this in the API.

    The ExperimentAxisQuery is a context manager and it is recommended that
    you use the following pattern::

        with ExperimentAxisQuery(...) as query:
            ...

    Lifecycle: maturing
    """

    def __init__(
        self,
        experiment: Experiment,
        measurement_name: str,
        *,
        obs_query: AxisQuery | None = None,
        var_query: AxisQuery | None = None,
        index_factory: IndexFactory = pd.Index,
    ) -> None:
        if measurement_name not in experiment.ms:
            raise ValueError("Measurement does not exist in the experiment")
        if obs_query is None:
            obs_query = AxisQuery()
        if var_query is None:
            var_query = AxisQuery()

        # Users often like to pass `foo=None` and we should let them
        obs_query = obs_query or AxisQuery()
        var_query = var_query or AxisQuery()

        self.experiment = experiment
        self.measurement_name = measurement_name

        self._matrix_axis_query = MatrixAxisQuery(obs=obs_query, var=var_query)
        self._joinids = JoinIDCache(self)
        self._indexer = AxisIndexer(
            self,
            index_factory=index_factory,
        )
        self._index_factory = index_factory

    def obs(
        self,
        *,
        column_names: Sequence[str] | None = None,
        batch_size: BatchSize | None = None,
        partitions: ReadPartitions | None = None,
        result_order: ResultOrderStr = _RO_AUTO,
        platform_config: PlatformConfig | None = None,
    ) -> ReadIter[pa.Table]:
        """Returns ``obs`` as an `Arrow table
        <https://arrow.apache.org/docs/python/generated/pyarrow.Table.html>`_
        iterator.

        Lifecycle: maturing
        """
        if batch_size is None:
            batch_size = BatchSize()
        obs_query = self._matrix_axis_query.obs
        return self._obs_df.read(
            obs_query.coords,
            value_filter=obs_query.value_filter,
            column_names=column_names,
            batch_size=batch_size,
            partitions=partitions,
            result_order=result_order,
            platform_config=platform_config,
        )

    def var(
        self,
        *,
        column_names: Sequence[str] | None = None,
        batch_size: BatchSize | None = None,
        partitions: ReadPartitions | None = None,
        result_order: ResultOrderStr = _RO_AUTO,
        platform_config: PlatformConfig | None = None,
    ) -> ReadIter[pa.Table]:
        """Returns ``var`` as an `Arrow table
        <https://arrow.apache.org/docs/python/generated/pyarrow.Table.html>`_
        iterator.

        Lifecycle: maturing
        """
        if batch_size is None:
            batch_size = BatchSize()
        var_query = self._matrix_axis_query.var
        return self._var_df.read(
            var_query.coords,
            value_filter=var_query.value_filter,
            column_names=column_names,
            batch_size=batch_size,
            partitions=partitions,
            result_order=result_order,
            platform_config=platform_config,
        )

    def obs_joinids(self) -> pa.IntegerArray:
        """Returns ``obs`` ``soma_joinids`` as an Arrow array.

        Lifecycle: maturing
        """
        return self._joinids.obs

    def var_joinids(self) -> pa.IntegerArray:
        """Returns ``var`` ``soma_joinids`` as an Arrow array.

        Lifecycle: maturing
        """
        return self._joinids.var

    @property
    def n_obs(self) -> int:
        """The number of ``obs`` axis query results.

        Lifecycle: maturing
        """
        return len(self.obs_joinids())

    @property
    def n_vars(self) -> int:
        """The number of ``var`` axis query results.

        Lifecycle: maturing
        """
        return len(self.var_joinids())

    @property
    def indexer(self) -> AxisIndexer:
        """A ``soma_joinid`` indexer for both ``obs`` and ``var`` axes.

        Lifecycle: maturing
        """
        return self._indexer

    def X(
        self,
        layer_name: str,
        *,
        batch_size: BatchSize | None = None,
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

        Lifecycle: maturing
        """
        if batch_size is None:
            batch_size = BatchSize()
        try:
            x_layer = self._ms.X[layer_name]
        except KeyError as ke:
            raise KeyError(f"{layer_name} is not present in X") from ke
        if not isinstance(x_layer, SparseNDArray):
            raise TypeError("X layers may only be sparse arrays")

        self._joinids.preload(self._threadpool)
        return x_layer.read(
            (self._joinids.obs, self._joinids.var),
            batch_size=batch_size,
            partitions=partitions,
            result_order=result_order,
            platform_config=platform_config,
        )

    def obsp(self, layer: str) -> SparseRead:
        """Returns an ``obsp`` layer as a sparse read.

        Lifecycle: maturing
        """
        joinids = self._joinids.obs
        return self._get_annotation_layer("obsp", layer).read((joinids, joinids))

    def varp(self, layer: str) -> SparseRead:
        """Returns a ``varp`` layer as a sparse read.

        Lifecycle: maturing
        """
        joinids = self._joinids.var
        return self._get_annotation_layer("varp", layer).read((joinids, joinids))

    def obsm(self, layer: str) -> SparseRead:
        """Returns an ``obsm`` layer as a sparse read.

        Lifecycle: maturing
        """
        return self._get_annotation_layer("obsm", layer).read((self._joinids.obs, slice(None)))

    def varm(self, layer: str) -> SparseRead:
        """Returns a ``varm`` layer as a sparse read.

        Lifecycle: maturing
        """
        return self._get_annotation_layer("varm", layer).read((self._joinids.var, slice(None)))

    def obs_scene_ids(self) -> pa.Array:
        """Returns a pyarrow array with scene ids that contain obs from this query.

        Lifecycle: experimental
        """
        try:
            obs_scene = self.experiment.obs_spatial_presence
        except KeyError as ke:
            raise KeyError("No obs_spatial_presence dataframe in this experiment.") from ke
        if not isinstance(obs_scene, DataFrame):
            raise TypeError(f"obs_spatial_presence must be a dataframe; got {type(obs_scene).__name__}.")

        full_table = obs_scene.read(
            coords=(self._joinids.obs, slice(None)),
            result_order=ResultOrder.COLUMN_MAJOR,
            value_filter="data != 0",
        ).concat()

        return pacomp.unique(full_table["scene_id"])

    def var_scene_ids(self) -> pa.Array:
        """Return a pyarrow array with scene ids that contain var from this
        query.

        Lifecycle: experimental
        """
        try:
            var_scene = self._ms.var_spatial_presence
        except KeyError as ke:
            raise KeyError(f"No var_spatial_presence dataframe in measurement '{self.measurement_name}'.") from ke
        if not isinstance(var_scene, DataFrame):
            raise TypeError(f"var_spatial_presence must be a dataframe; got {type(var_scene).__name__}.")

        full_table = var_scene.read(
            coords=(self._joinids.var, slice(None)),
            result_order=ResultOrder.COLUMN_MAJOR,
            value_filter="data != 0",
        ).concat()

        return pacomp.unique(full_table["scene_id"])

    def to_anndata(
        self,
        X_name: str | Sentinel | None = MISSING,
        *,
        column_names: AxisColumnNames | None = None,
        X_layers: Sequence[str] = (),
        obsm_layers: Sequence[str] = (),
        obsp_layers: Sequence[str] = (),
        varm_layers: Sequence[str] = (),
        varp_layers: Sequence[str] = (),
        drop_levels: bool = False,
        dask: SOMADaskConfig | None = None,
        obs_id_name: str | None = None,
        var_id_name: str | None = None,
    ) -> AnnData:
        """Exports the query to an in-memory ``AnnData`` object.

        Args:
            X_name:
                The X layer to read and return in the ``X`` slot.
                If unspecified (default), and the measurement contains an X layer named
                ``"data"``, it will be used. If ``None``, the returned AnnData will have
                ``X=None`` and ``layers`` will be unpopulated. If a string, that layer must
                exist in the measurement ``X`` collection.
            column_names:
                The columns in the ``var`` and ``obs`` dataframes to read.
            X_layers:
                Additional X layers to read and return in the ``layers`` slot.
            obsm_layers:
                Additional obsm layers to read and return in the obsm slot.
            obsp_layers:
                Additional obsp layers to read and return in the obsp slot.
            varm_layers:
                Additional varm layers to read and return in the varm slot.
            varp_layers:
                Additional varp layers to read and return in the varp slot.
            drop_levels:
                If true, drop unused categories from the ``obs`` and ``var`` dataframes.
                Defaults to ``False``.
            dask:
                If not ``None``, load the X layer as a Dask array. See
                :class:`DaskConfig` for details.
            obs_id_name:
                If specified, set this column as the index in the ``obs`` dataframe.
                If not specified, the default index column (see Notes below) will be determined,
                and set as the dataframe index.
            var_id_name:
                If specified, set this column as the index in the ``obs`` dataframe.
                If not specified, the default index column (see Notes below) will be determined,
                and set as the dataframe index.

        Notes:
        The default index column for the ``obs`` and ``var`` dataframes is determined with the following
        algorithm:

        - if the index column name is explicitly specified (e.g., via ``obs_id_name``), it
          will be used as the index.
        - if the original index column was written to the DataFrame metadata (e.g., during
          ingestion from AnnData), this will be used as the index.
        - if the dataframe contains a column named ``obs_id`` (for obs) or ``var_id`` (for var),
          this will be used as the index.
        - otherwise, the default Pandas index will be used.

        Lifecycle: experimental
        """
        if column_names is None:
            column_names = AxisColumnNames(obs=None, var=None)

        tp = self._threadpool
        x_collection = self._ms.X
        all_x_arrays: dict[str, SparseNDArray] = {}
        # Enforce identical semantics to io.to_anndata for X handling
        if X_name is None and X_layers:
            raise ValueError("If X_name is None, X_layers must not be provided")

        chosen_x: str | None = None
        if X_name is MISSING:
            if "data" in x_collection:
                chosen_x = "data"
        elif X_name is not None:
            if not isinstance(X_name, str) or not X_name:
                raise ValueError("X layer names must be specified as a string.")
            if X_name not in x_collection:
                raise ValueError(f"X layer name '{X_name}' not found in measurement")
            chosen_x = X_name

        # Gather X and any requested extra layers
        if chosen_x is not None:
            x_array = x_collection[chosen_x]
            if not isinstance(x_array, SparseNDArray):
                raise NotImplementedError("Dense array unsupported")
            all_x_arrays[chosen_x] = x_array
        for x_name in list(X_layers):
            if not isinstance(x_name, str) or not x_name:
                raise ValueError("X layer names must be specified as a string.")
            if x_name not in x_collection:
                raise ValueError("Unknown X layer name")
            x_array = x_collection[x_name]
            if not isinstance(x_array, SparseNDArray):
                raise NotImplementedError("Dense array unsupported")
            all_x_arrays[x_name] = x_array

        obs, var = tp.map(
            self._read_axis_dataframe,
            (AxisName.OBS, AxisName.VAR),
            (self._obs_df, self._var_df),
            (self._matrix_axis_query.obs, self._matrix_axis_query.var),
            (column_names, column_names),
            (obs_id_name, var_id_name),
            ("obs_id", "var_id"),
        )
        obs_joinids = self.obs_joinids()
        var_joinids = self.var_joinids()

        x_matrices = (
            {
                x_name: (
                    tp.submit(
                        _read_as_csr,
                        layer,
                        obs_joinids,
                        var_joinids,
                        self._indexer.by_obs,
                        self._indexer.by_var,
                    )
                    if not dask
                    else load_daskarray(
                        layer=layer,
                        coords=(obs_joinids, var_joinids),
                        **dask,
                    )
                )
                for x_name, layer in all_x_arrays.items()
            }
            if (X_name is not None and X_name is not MISSING) or chosen_x is not None or X_layers
            else {}
        )
        x_future = x_matrices.pop(chosen_x) if chosen_x is not None else None

        obsm_future = {
            key: tp.submit(
                _read_inner_ndarray,
                self._get_annotation_layer("obsm", key),
                obs_joinids,
                self.indexer.by_obs,
            )
            for key in obsm_layers
        }
        varm_future = {
            key: tp.submit(
                _read_inner_ndarray,
                self._get_annotation_layer("varm", key),
                var_joinids,
                self.indexer.by_var,
            )
            for key in varm_layers
        }
        obsp_future = {
            key: tp.submit(
                _read_as_csr,
                self._get_annotation_layer("obsp", key),
                obs_joinids,
                obs_joinids,
                self.indexer.by_obs,
                self.indexer.by_obs,
            )
            for key in obsp_layers
        }
        varp_future = {
            key: tp.submit(
                _read_as_csr,
                self._get_annotation_layer("varp", key),
                var_joinids,
                var_joinids,
                self.indexer.by_var,
                self.indexer.by_var,
            )
            for key in varp_layers
        }

        # Drop unused categories on axis dataframes if requested
        if drop_levels:
            for name in obs:
                if obs[name].dtype.name == "category":
                    obs[name] = obs[name].cat.remove_unused_categories()
            for name in var:
                if var[name].dtype.name == "category":
                    var[name] = var[name].cat.remove_unused_categories()

        return AnnData(
            X=(x_future.result() if isinstance(x_future, Future) else x_future) if chosen_x is not None else None,
            obs=obs,
            var=var,
            obsm=(_resolve_futures(obsm_future) or None),
            obsp=(_resolve_futures(obsp_future) or None),
            varm=(_resolve_futures(varm_future) or None),
            varp=(_resolve_futures(varp_future) or None),
            layers=_resolve_futures(x_matrices),
        )

    def to_spatialdata(  # type: ignore[no-untyped-def]  # noqa: ANN202
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
        scene_presence_mode: str = "obs",
        dask: SOMADaskConfig | None = None,
    ):
        """Returns a SpatialData object containing the query results.

        This is a low-level routine intended to be used by loaders for other
        in-core formats, such as AnnData, which can be created from the
        resulting objects.

        Args:
            X_name: The X layer to read and return in the ``X`` slot.
            column_names: The columns in the ``var`` and ``obs`` dataframes
                to read.
            X_layers: Additional X layers to read and return in the ``layers`` slot.
            obsm_layers: Additional obsm layers to read and return in the obsm slot.
            obsp_layers: Additional obsp layers to read and return in the obsp slot.
            varm_layers: Additional varm layers to read and return in the varm slot.
            varp_layers: Additional varp layers to read and return in the varp slot.
            drop_levels: If ``True`` remove unused categories from the AnnData table
                with the measurement data.
            scene_presence_mode: Method for determining what scenes to return data
                from. Valid options are ``obs`` (use ``obs_spatial_presence``
                dataframe) and ``var`` (use ``var_spatial_presence`` dataframe).
                Defaults to ``obs``.
        """
        from .io.spatial._spatialdata_util import _spatial_to_spatialdata

        warnings.warn(SPATIAL_DISCLAIMER, stacklevel=2)

        # Get a list of scenes to add to SpatialData object.
        if scene_presence_mode == "obs":
            scene_names = tuple(str(scene_name) for scene_name in self.obs_scene_ids())
        elif scene_presence_mode == "var":
            scene_names = tuple(str(scene_name) for scene_name in self.var_scene_ids())
        else:
            raise ValueError(f"Invalid scene presence mode '{scene_presence_mode}'. Valid options are 'obs' and 'var'.")

        # Get the anndata table.
        ad = self.to_anndata(
            X_name,
            column_names=column_names,
            X_layers=X_layers,
            obsm_layers=obsm_layers,
            obsp_layers=obsp_layers,
            varm_layers=varm_layers,
            varp_layers=varp_layers,
            drop_levels=drop_levels,
            dask=dask,
        )

        return _spatial_to_spatialdata(
            self.experiment.spatial,
            {self.measurement_name: ad},
            scene_names=scene_names,
            measurement_names=(self.measurement_name,),
        )

    # Context management.
    # Currently a no-op, but, part of public API so we retain it.
    def close(self) -> None:
        pass

    def __enter__(self) -> Self:
        return self

    def __exit__(self, *_: Any) -> None:  # noqa: ANN401
        pass

    # Internals

    def _read_axis_dataframe(
        self,
        axis: AxisName,
        axis_df: DataFrame,
        axis_query: AxisQuery,
        axis_column_names: AxisColumnNames,
        default_index_name: str | None = None,
        fallback_index_name: str | None = None,
    ) -> pd.DataFrame:
        """Reads the specified axis. Will cache join IDs if not present."""
        column_names = axis_column_names.get(axis.value)

        # If we can cache join IDs, prepare to add them to the cache.
        joinids_cached = self._joinids._is_cached(axis)
        query_columns = column_names
        added_soma_joinid_to_columns = False
        if not joinids_cached and column_names is not None and "soma_joinid" not in column_names:
            # If we want to fill the join ID cache, ensure that we query the
            # soma_joinid column so that it is included in the results.
            # We'll filter it out later.
            query_columns = ["soma_joinid", *list(column_names)]
            added_soma_joinid_to_columns = True

        # Do the actual query.
        arrow_table = axis_df.read(
            coords=axis_query.coords,
            value_filter=axis_query.value_filter,
            column_names=query_columns,
        ).concat()

        # Update the cache if needed. We can do this because no matter what
        # other columns are queried for, the contents of the ``soma_joinid``
        # column will be the same and can be safely stored.
        if not joinids_cached:
            setattr(
                self._joinids,
                axis.value,
                arrow_table.column("soma_joinid").combine_chunks(),
            )

        # Drop soma_joinid column if we added it solely for use in filling
        # the joinid cache.
        if added_soma_joinid_to_columns:
            arrow_table = arrow_table.drop(["soma_joinid"])

        # Read and validate the "original index metadata" stored alongside this SOMA DataFrame.
        original_index_metadata = json.loads(axis_df.metadata.get(SOMA_DATAFRAME_ORIGINAL_INDEX_NAME_JSON, "null"))
        if not (original_index_metadata is None or isinstance(original_index_metadata, str)):
            raise ValueError(
                f"{axis_df.uri}: invalid {SOMA_DATAFRAME_ORIGINAL_INDEX_NAME_JSON} metadata: {original_index_metadata}",
            )

        df: pd.DataFrame = arrow_table.to_pandas()

        default_index_name = default_index_name or original_index_metadata
        _df_set_index(df, default_index_name, fallback_index_name)

        return df

    def _get_annotation_layer(self, annotation_name: str, layer_name: str) -> SparseNDArray:
        """Helper function to make error messages consistent.

        Args:
            annotation_name:
                Name of the annotation (e.g. obsm, varp).
            layer_name:
                Name of the layer.
        """
        try:
            coll = self._ms[annotation_name]
        except KeyError:
            raise ValueError(f"Measurement does not contain {annotation_name!r} data.") from None
        if not isinstance(coll, Collection):
            raise TypeError(f"Unexpected SOMA type {type(coll).__name__} for {annotation_name!r}.")

        try:
            layer = coll[layer_name]
        except KeyError:
            raise ValueError(f"layer {layer_name!r} is not available in {annotation_name!r}.") from None
        if not isinstance(layer, SparseNDArray):
            raise TypeError(
                f"Unexpected SOMA type {type(layer).__name__} stored in {annotation_name!r} layer {layer_name!r}.",
            )
        return layer

    @property
    def _obs_df(self) -> DataFrame:
        return self.experiment.obs

    @property
    def _ms(self) -> Measurement:
        return self.experiment.ms[self.measurement_name]

    @property
    def _var_df(self) -> DataFrame:
        return self._ms.var

    @property
    def _threadpool(self) -> ThreadPoolExecutor:
        """Returns the threadpool provided by the experiment's context.
        If not available, creates a thread pool just in time.
        """
        return self.experiment.context.threadpool


@attrs.define(frozen=True)
class MatrixAxisQuery:
    """The per-axis user query definition. Private."""

    obs: AxisQuery
    var: AxisQuery


@attrs.define
class JoinIDCache:
    """A cache for per-axis join ids in the query. Private."""

    owner: ExperimentAxisQuery

    _cached_obs: pa.IntegerArray | None = None
    _cached_var: pa.IntegerArray | None = None

    def _is_cached(self, axis: AxisName) -> bool:
        field = "_cached_" + axis.value
        return getattr(self, field) is not None

    def preload(self, pool: ThreadPoolExecutor) -> None:
        if self._cached_obs is not None and self._cached_var is not None:
            return
        obs_ft = pool.submit(lambda: self.obs)
        var_ft = pool.submit(lambda: self.var)
        # Wait for them and raise in case of error.
        obs_ft.result()
        var_ft.result()

    @property
    def obs(self) -> pa.IntegerArray:
        """Join IDs for the obs axis. Will load and cache if not already."""
        if not self._cached_obs:
            self._cached_obs = load_joinids(self.owner._obs_df, self.owner._matrix_axis_query.obs)
        return self._cached_obs

    @obs.setter
    def obs(self, val: pa.IntegerArray) -> None:
        self._cached_obs = val

    @property
    def var(self) -> pa.IntegerArray:
        """Join IDs for the var axis. Will load and cache if not already."""
        if not self._cached_var:
            self._cached_var = load_joinids(self.owner._var_df, self.owner._matrix_axis_query.var)
        return self._cached_var

    @var.setter
    def var(self, val: pa.IntegerArray) -> None:
        self._cached_var = val


def load_joinids(df: DataFrame, axq: AxisQuery) -> pa.IntegerArray:
    tbl = df.read(
        axq.coords,
        value_filter=axq.value_filter,
        column_names=["soma_joinid"],
    ).concat()
    return tbl.column("soma_joinid").combine_chunks()


def _read_inner_ndarray(
    matrix: SparseNDArray,
    joinids: pa.IntegerArray,
    indexer: Callable[[Numpyable], npt.NDArray[np.intp]],
) -> npt.NDArray[np.float32]:
    table = matrix.read((joinids, slice(None))).tables().concat()
    dtype = matrix.schema.field("soma_data").type.to_pandas_dtype()

    idx = indexer(table["soma_dim_0"])
    n_row = len(joinids)

    # Prior to the shape / current_domain change in release 1.15.0, there
    # was no way to determine the "max" number of features, aka shape, of
    # the sparse obsm/varm array. Guess based upon the array contents.
    n_col = matrix.shape[1] if matrix.tiledbsoma_has_upgraded_shape else pa.compute.max(table["soma_dim_1"]).as_py() + 1

    z: npt.NDArray[np.float32] = np.zeros(n_row * n_col, dtype=dtype)
    np.put(z, idx * n_col + table["soma_dim_1"], table["soma_data"])
    return z.reshape(n_row, n_col)


def _read_as_csr(
    matrix: SparseNDArray,
    d0_joinids_arr: pa.IntegerArray,
    d1_joinids_arr: pa.IntegerArray,
    d0_indexer: Callable[[Numpyable], npt.NDArray[np.intp]],
    d1_indexer: Callable[[Numpyable], npt.NDArray[np.intp]],
) -> sp.csr_matrix:
    d0_joinids = d0_joinids_arr.to_numpy()
    d1_joinids = d1_joinids_arr.to_numpy()
    try:
        # frag_cell_count is >= nnz, as it does not account for deletes and double-counts updates
        frag_cell_count: int | None = matrix._handle.fragment_cell_count()
    except SOMAError:
        frag_cell_count = None

    # if able, downcast from int64 - reduces working memory
    index_dtype = np.int32 if max(len(d0_joinids), len(d1_joinids)) < np.iinfo(np.int32).max else np.int64
    pa_schema = pa.schema(
        [
            pa.field("soma_dim_0", pa.from_numpy_dtype(index_dtype)),
            pa.field("soma_dim_1", pa.from_numpy_dtype(index_dtype)),
            matrix.schema.field("soma_data"),
        ],
    )

    def _read_and_reindex(X: SparseNDArray, oids: npt.NDArray[np.int64], vids: npt.NDArray[np.int64]) -> pa.Table:
        def _reindex(batch: pa.RecordBatch) -> pa.RecordBatch:
            return pa.RecordBatch.from_pydict(
                {
                    "soma_dim_0": d0_indexer(batch["soma_dim_0"]).astype(index_dtype),
                    "soma_dim_1": d1_indexer(batch["soma_dim_1"]).astype(index_dtype),
                    "soma_data": batch["soma_data"],
                },
                schema=pa_schema,
            )

        return pa.Table.from_batches(
            (_reindex(batch) for tbl in X.read(coords=(oids, vids)).tables() for batch in tbl.to_batches()),
            schema=pa_schema,
        )

    approx_X_shape = tuple(b - a + 1 for a, b in matrix.non_empty_domain())
    # heuristically derived number (benchmarking). Thesis is that this is roughly 80% of a 1 GiB io buffer,
    # which is the default for SOMA. If we have fast NNZ, use to partition based upon memory size. If we do not
    # have fast NNZ, pick a default partition size which is large, but not so large as to rule out reasonable
    # parallelism (in this case, calculated based on typical scRNASeq assay density of 3-6K features).
    target_point_count = 96 * 1024**2
    fallback_row_count = 32768
    # compute partition size from array density and target point count, rounding to nearest 1024.
    partition_size = (
        max(
            1024 * round(approx_X_shape[0] * target_point_count / frag_cell_count / 1024),
            1024,
        )
        if frag_cell_count is not None and frag_cell_count > 0
        else min(fallback_row_count, approx_X_shape[0])
    )
    splits = list(
        range(
            partition_size,
            len(d0_joinids) - partition_size + 1,
            partition_size,
        ),
    )
    if len(splits) > 1:
        d0_joinids_splits = np.array_split(np.partition(d0_joinids, splits), splits)
        tp = matrix.context.threadpool
        tbl = pa.concat_tables(
            tp.map(
                _read_and_reindex,
                (matrix,) * len(d0_joinids_splits),
                d0_joinids_splits,
                (d1_joinids,) * len(d0_joinids_splits),
            ),
        )

    else:
        tbl = _read_and_reindex(matrix, d0_joinids, d1_joinids)

    return CompressedMatrix.from_soma(tbl, (len(d0_joinids), len(d1_joinids)), "csr", True, matrix.context).to_scipy()
