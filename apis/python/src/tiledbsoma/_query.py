# Copyright (c) 2021-2023 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2023 TileDB, Inc.
#
# Licensed under the MIT License.

"""Implementation of a SOMA Experiment.
"""

from __future__ import annotations

import enum
import warnings
from concurrent.futures import ThreadPoolExecutor
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Dict,
    Literal,
    Mapping,
    Protocol,
    Sequence,
    TypeVar,
    cast,
    overload,
)

import attrs
import numpy as np
import numpy.typing as npt
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pacomp
import scipy.sparse as sp
from anndata import AnnData
from somacore import (
    AxisQuery,
    DataFrame,
    ReadIter,
    SparseRead,
    query,
)
from somacore.data import _RO_AUTO
from somacore.options import (
    BatchSize,
    PlatformConfig,
    ReadPartitions,
    ResultOrder,
    ResultOrderStr,
)
from somacore.query.query import (
    AxisColumnNames,
    Numpyable,
)
from somacore.query.types import IndexFactory, IndexLike
from typing_extensions import Self

if TYPE_CHECKING:
    from ._experiment import Experiment
from ._constants import SPATIAL_DISCLAIMER
from ._fastercsx import CompressedMatrix
from ._measurement import Measurement
from ._sparse_nd_array import SparseNDArray
from ._util import _resolve_futures

_T = TypeVar("_T")
_T_co = TypeVar("_T_co", covariant=True)


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
        return super().value  # type: ignore[no-any-return]

    @overload
    def getattr_from(self, __source: _HasObsVar[_T]) -> _T: ...

    @overload
    def getattr_from(
        self, __source: Any, *, pre: Literal[""], suf: Literal[""]
    ) -> object: ...

    @overload
    def getattr_from(
        self, __source: Any, *, pre: str = ..., suf: str = ...
    ) -> object: ...

    def getattr_from(self, __source: Any, *, pre: str = "", suf: str = "") -> object:
        """Equivalent to ``something.<pre><obs/var><suf>``."""
        return getattr(__source, pre + self.value + suf)

    def getitem_from(
        self, __source: Mapping[str, "_T"], *, pre: str = "", suf: str = ""
    ) -> _T:
        """Equivalent to ``something[pre + "obs"/"var" + suf]``."""
        return __source[pre + self.value + suf]


@attrs.define
class AxisIndexer(query.AxisIndexer):
    """
    Given a query, provides index-building services for obs/var axis.

    Lifecycle: maturing
    """

    query: "ExperimentAxisQuery"
    _index_factory: IndexFactory
    _cached_obs: IndexLike | None = None
    _cached_var: IndexLike | None = None

    @property
    def _obs_index(self) -> IndexLike:
        """Private. Return an index for the ``obs`` axis."""
        if self._cached_obs is None:
            self._cached_obs = self._index_factory(self.query.obs_joinids().to_numpy())
        return self._cached_obs

    @property
    def _var_index(self) -> IndexLike:
        """Private. Return an index for the ``var`` axis."""
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


@attrs.define(frozen=True)
class AxisQueryResult:
    """The result of running :meth:`ExperimentAxisQuery.read`. Private."""

    obs: pd.DataFrame
    """Experiment.obs query slice, as a pandas DataFrame"""
    var: pd.DataFrame
    """Experiment.ms[...].var query slice, as a pandas DataFrame"""
    X: sp.csr_matrix
    """Experiment.ms[...].X[...] query slice, as a SciPy sparse.csr_matrix """
    X_layers: Dict[str, sp.csr_matrix] = attrs.field(factory=dict)
    """Any additional X layers requested, as SciPy sparse.csr_matrix(s)"""
    obsm: Dict[str, npt.NDArray[Any]] = attrs.field(factory=dict)
    """Experiment.obsm query slice, as a numpy ndarray"""
    obsp: Dict[str, sp.csr_matrix] = attrs.field(factory=dict)
    """Experiment.obsp query slice, as SciPy sparse.csr_matrix(s)"""
    varm: Dict[str, npt.NDArray[Any]] = attrs.field(factory=dict)
    """Experiment.varm query slice, as a numpy ndarray"""
    varp: Dict[str, sp.csr_matrix] = attrs.field(factory=dict)
    """Experiment.varp query slice, as SciPy sparse.csr_matrix(s)"""

    def to_anndata(self) -> AnnData:
        return AnnData(
            X=self.X,
            obs=self.obs,
            var=self.var,
            obsm=(self.obsm or None),
            obsp=(self.obsp or None),
            varm=(self.varm or None),
            varp=(self.varp or None),
            layers=(self.X_layers or None),
        )


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
        experiment: "Experiment",
        measurement_name: str,
        *,
        obs_query: AxisQuery = AxisQuery(),
        var_query: AxisQuery = AxisQuery(),
        index_factory: IndexFactory = pd.Index,
    ):
        if measurement_name not in experiment.ms:
            raise ValueError("Measurement does not exist in the experiment")

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
        batch_size: BatchSize = BatchSize(),
        partitions: ReadPartitions | None = None,
        result_order: ResultOrderStr = _RO_AUTO,
        platform_config: PlatformConfig | None = None,
    ) -> ReadIter[pa.Table]:
        """Returns ``obs`` as an `Arrow table
        <https://arrow.apache.org/docs/python/generated/pyarrow.Table.html>`_
        iterator.

        Lifecycle: maturing
        """
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
        batch_size: BatchSize = BatchSize(),
        partitions: ReadPartitions | None = None,
        result_order: ResultOrderStr = _RO_AUTO,
        platform_config: PlatformConfig | None = None,
    ) -> ReadIter[pa.Table]:
        """Returns ``var`` as an `Arrow table
        <https://arrow.apache.org/docs/python/generated/pyarrow.Table.html>`_
        iterator.

        Lifecycle: maturing
        """
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
        batch_size: BatchSize = BatchSize(),
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
        return self._axisp_get_array(AxisName.OBS, layer).read((joinids, joinids))

    def varp(self, layer: str) -> SparseRead:
        """Returns a ``varp`` layer as a sparse read.

        Lifecycle: maturing
        """
        joinids = self._joinids.var
        return self._axisp_get_array(AxisName.VAR, layer).read((joinids, joinids))

    def obsm(self, layer: str) -> SparseRead:
        """Returns an ``obsm`` layer as a sparse read.
        Lifecycle: maturing
        """
        return self._axism_get_array(AxisName.OBS, layer).read(
            (self._joinids.obs, slice(None))
        )

    def varm(self, layer: str) -> SparseRead:
        """Returns a ``varm`` layer as a sparse read.
        Lifecycle: maturing
        """
        return self._axism_get_array(AxisName.VAR, layer).read(
            (self._joinids.var, slice(None))
        )

    def obs_scene_ids(self) -> pa.Array:
        """Returns a pyarrow array with scene ids that contain obs from this
        query.

        Lifecycle: experimental
        """
        try:
            obs_scene = self.experiment.obs_spatial_presence
        except KeyError as ke:
            raise KeyError(
                "No obs_spatial_presence dataframe in this experiment."
            ) from ke
        if not isinstance(obs_scene, DataFrame):
            raise TypeError(
                f"obs_spatial_presence must be a dataframe; got "
                f"{type(obs_scene).__name__}."
            )

        full_table = obs_scene.read(
            coords=((AxisName.OBS.getattr_from(self._joinids), slice(None))),
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
            raise KeyError(
                f"No var_spatial_presence dataframe in measurement "
                f"'{self.measurement_name}'."
            ) from ke
        if not isinstance(var_scene, DataFrame):
            raise TypeError(
                f"var_spatial_presence must be a dataframe; got "
                f"{type(var_scene).__name__}."
            )

        full_table = var_scene.read(
            coords=((AxisName.VAR.getattr_from(self._joinids), slice(None))),
            result_order=ResultOrder.COLUMN_MAJOR,
            value_filter="data != 0",
        ).concat()

        return pacomp.unique(full_table["scene_id"])

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
        ad = self._read(
            X_name,
            column_names=column_names or AxisColumnNames(obs=None, var=None),
            X_layers=X_layers,
            obsm_layers=obsm_layers,
            obsp_layers=obsp_layers,
            varm_layers=varm_layers,
            varp_layers=varp_layers,
        ).to_anndata()

        # Drop unused categories on axis dataframes if requested
        if drop_levels:
            for name in ad.obs:
                if ad.obs[name].dtype.name == "category":
                    ad.obs[name] = ad.obs[name].cat.remove_unused_categories()
            for name in ad.var:
                if ad.var[name].dtype.name == "category":
                    ad.var[name] = ad.var[name].cat.remove_unused_categories()

        return ad

    def to_spatialdata(  # type: ignore[no-untyped-def]
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
    ):
        """Returns a SpatialData object containing the query results

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

        from spatialdata import SpatialData

        from .io.spatial.outgest import _add_scene_to_spatialdata

        warnings.warn(SPATIAL_DISCLAIMER)

        # Get a list of scenes to add to SpatialData object.
        if scene_presence_mode == "obs":
            scene_ids = self.obs_scene_ids()
        elif scene_presence_mode == "var":
            scene_ids = self.var_scene_ids()
        else:
            raise ValueError(
                f"Invalid scene presence mode '{scene_presence_mode}'. Valid options "
                f"are 'obs' and 'var'."
            )

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
        )
        sdata = SpatialData(tables={self.measurement_name: ad})

        for scene_id in scene_ids:
            scene = self.experiment.spatial[str(scene_id)]
            _add_scene_to_spatialdata(
                sdata,
                scene_id=str(scene_id),
                scene=scene,
                obs_id_name="soma_joinid",
                var_id_name="soma_joinid",
                measurement_names=(self.measurement_name,),
            )

        return sdata

    # Context management

    def __enter__(self) -> Self:
        return self

    def __exit__(self, *_: Any) -> None:
        pass

    # Internals

    def _read(
        self,
        X_name: str,
        *,
        column_names: AxisColumnNames,
        X_layers: Sequence[str],
        obsm_layers: Sequence[str] = (),
        obsp_layers: Sequence[str] = (),
        varm_layers: Sequence[str] = (),
        varp_layers: Sequence[str] = (),
    ) -> AxisQueryResult:
        """Reads the entire query result in memory.

        This is a low-level routine intended to be used by loaders for other
        in-core formats, such as AnnData, which can be created from the
        resulting objects.

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
        """
        tp = self._threadpool
        x_collection = self._ms.X
        all_x_names = [X_name] + list(X_layers)
        all_x_arrays: Dict[str, SparseNDArray] = {}
        for _xname in all_x_names:
            if not isinstance(_xname, str) or not _xname:
                raise ValueError("X layer names must be specified as a string.")
            if _xname not in x_collection:
                raise ValueError("Unknown X layer name")
            x_array = x_collection[_xname]
            if not isinstance(x_array, SparseNDArray):
                raise NotImplementedError("Dense array unsupported")
            all_x_arrays[_xname] = x_array

        obs_table, var_table = tp.map(
            self._read_axis_dataframe,
            (AxisName.OBS, AxisName.VAR),
            (column_names, column_names),
        )
        obs_joinids = self.obs_joinids()
        var_joinids = self.var_joinids()

        x_matrices = {
            _xname: tp.submit(
                _read_as_csr,
                layer,
                obs_joinids,
                var_joinids,
                self._indexer.by_obs,
                self._indexer.by_var,
            )
            for _xname, layer in all_x_arrays.items()
        }
        x_future = x_matrices.pop(X_name)

        obsm_future = {
            key: tp.submit(self._axism_inner_ndarray, AxisName.OBS, key)
            for key in obsm_layers
        }
        varm_future = {
            key: tp.submit(self._axism_inner_ndarray, AxisName.VAR, key)
            for key in varm_layers
        }
        obsp_future = {
            key: tp.submit(self._axisp_inner_sparray, AxisName.OBS, key)
            for key in obsp_layers
        }
        varp_future = {
            key: tp.submit(self._axisp_inner_sparray, AxisName.VAR, key)
            for key in varp_layers
        }

        obs = obs_table.to_pandas()
        obs.index = obs.index.astype(str)

        var = var_table.to_pandas()
        var.index = var.index.astype(str)

        return AxisQueryResult(
            obs=obs,
            var=var,
            X=x_future.result(),
            obsm=_resolve_futures(obsm_future),
            obsp=_resolve_futures(obsp_future),
            varm=_resolve_futures(varm_future),
            varp=_resolve_futures(varp_future),
            X_layers=_resolve_futures(x_matrices),
        )

    def _read_axis_dataframe(
        self,
        axis: AxisName,
        axis_column_names: AxisColumnNames,
    ) -> pa.Table:
        """Reads the specified axis. Will cache join IDs if not present."""
        column_names = axis_column_names.get(axis.value)

        axis_df = axis.getattr_from(self, pre="_", suf="_df")
        assert isinstance(axis_df, DataFrame)
        axis_query = axis.getattr_from(self._matrix_axis_query)

        # If we can cache join IDs, prepare to add them to the cache.
        joinids_cached = self._joinids._is_cached(axis)
        query_columns = column_names
        added_soma_joinid_to_columns = False
        if (
            not joinids_cached
            and column_names is not None
            and "soma_joinid" not in column_names
        ):
            # If we want to fill the join ID cache, ensure that we query the
            # soma_joinid column so that it is included in the results.
            # We'll filter it out later.
            query_columns = ["soma_joinid"] + list(column_names)
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
        return arrow_table

    def _axisp_get_array(
        self,
        axis: AxisName,
        layer: str,
    ) -> SparseNDArray:
        p_name = f"{axis.value}p"
        try:
            ms = self._ms
            axisp = ms.obsp if axis.value == "obs" else ms.varp
        except (AttributeError, KeyError):
            raise ValueError(f"Measurement does not contain {p_name} data")

        try:
            axisp_layer = axisp[layer]
        except KeyError:
            raise ValueError(f"layer {layer!r} is not available in {p_name}")
        if not isinstance(axisp_layer, SparseNDArray):
            raise TypeError(
                f"Unexpected SOMA type {type(axisp_layer).__name__}"
                f" stored in {p_name} layer {layer!r}"
            )

        return axisp_layer

    def _axism_get_array(
        self,
        axis: AxisName,
        layer: str,
    ) -> SparseNDArray:
        m_name = f"{axis.value}m"

        try:
            ms = self._ms
            axism = ms.obsm if axis.value == "obs" else ms.varm
        except (AttributeError, KeyError):
            raise ValueError(f"Measurement does not contain {m_name} data")

        try:
            axism_layer = axism[layer]
        except KeyError:
            raise ValueError(f"layer {layer!r} is not available in {m_name}")

        if not isinstance(axism_layer, SparseNDArray):
            raise TypeError(f"Unexpected SOMA type stored in '{m_name}' layer")

        return axism_layer

    def _convert_to_ndarray(
        self, axis: AxisName, table: pa.Table, n_row: int, n_col: int
    ) -> npt.NDArray[np.float32]:
        indexer = cast(
            Callable[[Numpyable], npt.NDArray[np.intp]],
            axis.getattr_from(self.indexer, pre="by_"),
        )
        idx = indexer(table["soma_dim_0"])
        z: npt.NDArray[np.float32] = np.zeros(n_row * n_col, dtype=np.float32)
        np.put(z, idx * n_col + table["soma_dim_1"], table["soma_data"])
        return z.reshape(n_row, n_col)

    def _axisp_inner_sparray(
        self,
        axis: AxisName,
        layer: str,
    ) -> sp.csr_matrix:
        joinids = axis.getattr_from(self._joinids)
        indexer = cast(
            Callable[[Numpyable], npt.NDArray[np.intp]],
            axis.getattr_from(self.indexer, pre="by_"),
        )
        return _read_as_csr(
            self._axisp_get_array(axis, layer), joinids, joinids, indexer, indexer
        )

    def _axism_inner_ndarray(
        self,
        axis: AxisName,
        layer: str,
    ) -> npt.NDArray[np.float32]:
        joinids = axis.getattr_from(self._joinids)
        table = (
            self._axism_get_array(axis, layer)
            .read((joinids, slice(None)))
            .tables()
            .concat()
        )

        n_row = len(joinids)
        n_col = len(table["soma_dim_1"].unique())

        return self._convert_to_ndarray(axis, table, n_row, n_col)

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
        """
        Returns the threadpool provided by the experiment's context.
        If not available, creates a thread pool just in time."""
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
            self._cached_obs = load_joinids(
                self.owner._obs_df, self.owner._matrix_axis_query.obs
            )
        return self._cached_obs

    @obs.setter
    def obs(self, val: pa.IntegerArray) -> None:
        self._cached_obs = val

    @property
    def var(self) -> pa.IntegerArray:
        """Join IDs for the var axis. Will load and cache if not already."""
        if not self._cached_var:
            self._cached_var = load_joinids(
                self.owner._var_df, self.owner._matrix_axis_query.var
            )
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


def _read_as_csr(
    matrix: SparseNDArray,
    d0_joinids_arr: pa.IntegerArray,
    d1_joinids_arr: pa.IntegerArray,
    d0_indexer: Callable[[Numpyable], npt.NDArray[np.intp]],
    d1_indexer: Callable[[Numpyable], npt.NDArray[np.intp]],
) -> sp.csr_matrix:

    d0_joinids = d0_joinids_arr.to_numpy()
    d1_joinids = d1_joinids_arr.to_numpy()
    nnz = matrix.nnz

    # if able, downcast from int64 - reduces working memory
    index_dtype = (
        np.int32
        if max(len(d0_joinids), len(d1_joinids)) < np.iinfo(np.int32).max
        else np.int64
    )
    pa_schema = pa.schema(
        [
            pa.field("soma_dim_0", pa.from_numpy_dtype(index_dtype)),
            pa.field("soma_dim_1", pa.from_numpy_dtype(index_dtype)),
            matrix.schema.field("soma_data"),
        ]
    )

    def _read_and_reindex(
        X: SparseNDArray, oids: npt.NDArray[np.int64], vids: npt.NDArray[np.int64]
    ) -> pa.Table:
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
            (
                _reindex(_batch)
                for _tbl in X.read(coords=(oids, vids)).tables()
                for _batch in _tbl.to_batches()
            ),
            schema=pa_schema,
        )

    approx_X_shape = tuple(b - a + 1 for a, b in matrix.non_empty_domain())
    # heuristically derived number (benchmarking). Thesis is that this is roughly 80% of a 1 GiB io buffer,
    # which is the default for SOMA.
    target_point_count = 96 * 1024**2
    # compute partition size from array density and target point count, rounding to nearest 1024.
    partition_size = (
        max(1024 * round(approx_X_shape[0] * target_point_count / nnz / 1024), 1024)
        if nnz > 0
        else approx_X_shape[0]
    )
    splits = list(
        range(
            partition_size,
            len(d0_joinids) - partition_size + 1,
            partition_size,
        )
    )
    if len(splits) > 0:
        d0_joinids_splits = np.array_split(np.partition(d0_joinids, splits), splits)
        tp = matrix.context.threadpool
        tbl = pa.concat_tables(
            tp.map(
                _read_and_reindex,
                (matrix,) * len(d0_joinids_splits),
                d0_joinids_splits,
                (d1_joinids,) * len(d0_joinids_splits),
            )
        )

    else:
        tbl = _read_and_reindex(matrix, d0_joinids, d1_joinids)

    return CompressedMatrix.from_soma(
        tbl, (len(d0_joinids), len(d1_joinids)), "csr", True, matrix.context
    ).to_scipy()
