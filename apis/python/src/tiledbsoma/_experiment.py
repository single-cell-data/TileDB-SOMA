# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""Implementation of a SOMA Experiment."""

from __future__ import annotations

import functools
from typing import ClassVar

import attrs
import pyarrow as pa
from somacore import experiment, options, query

from . import pytiledbsoma as clib
from ._collection import Collection, CollectionBase
from ._dataframe import DataFrame
from ._dense_nd_array import DenseNDArray
from ._exception import SOMAError
from ._geometry_dataframe import GeometryDataFrame
from ._indexer import IntIndexer
from ._measurement import Measurement
from ._point_cloud_dataframe import PointCloudDataFrame
from ._query import ExperimentAxisQuery
from ._scene import Scene
from ._soma_context import SOMAContext
from ._soma_object import SOMAObject
from ._sparse_nd_array import SparseNDArray
from .options._tiledb_create_write_options import TileDBDeleteOptions


class Experiment(
    CollectionBase[SOMAObject],
    experiment.Experiment[
        DataFrame,
        Collection[Measurement],  # type: ignore[type-var]
        Collection[Scene],  # type: ignore[type-var]
        SOMAObject,
    ],
):
    """A collection subtype that combines observations and measurements
    from an individual experiment.

    In single cell biology, this can represent multiple modes of measurement
    across a single collection of cells (i.e., a "multimodal dataset").
    Within an experiment, a set of measurements on a single set of variables
    (i.e., features) is represented as a :class:`Measurement`.

    Attributes:
        obs (DataFrame):
            Primary annotations on the observation axis. The contents of the
            ``soma_joinid`` column define the observation index domain,
            AKA ``obs_id``. All observations for the Experiment must be
            defined in this dataframe.
        ms (Collection):
            A collection of named measurements.
        spatial (Collection):
            A collection of spatial scenes.

    Example:
        >>> import tiledbsoma
        >>> with tiledbsoma.open("/path/to/experiment") as exp:
        ...     # While users can interact directly with an Experiment's fields:
        ...     obs_df = exp.obs
        ...
        ...     # the primary use case is to run queries on the experiment data.
        ...     q = exp.axis_query(
        ...         "mtdna",
        ...         obs_query=tiledbsoma.AxisQuery(value_filter="tissue == 'lung'"),
        ...         var_query=tiledbsoma.AxisQuery(coords=(slice(50, 100),)),
        ...     )
        ...     query_obs = q.obs().concat().to_pandas()
        ...     query_var = q.var().concat().to_pandas()

    Lifecycle:
        Maturing.
    """

    __slots__ = ()
    _handle_type = clib.SOMAExperiment

    _subclass_constrained_soma_types: ClassVar[dict[str, tuple[str, ...]]] = {
        "obs": ("SOMADataFrame",),
        "ms": ("SOMACollection",),
        "spatial": ("SOMACollection",),
        "obs_spatial_presence": ("SOMADataFrame",),
    }

    def axis_query(
        self,
        measurement_name: str,
        *,
        obs_query: query.AxisQuery | None = None,
        var_query: query.AxisQuery | None = None,
    ) -> ExperimentAxisQuery:
        """Creates an axis query over this experiment.
        Lifecycle: Maturing.
        """
        return ExperimentAxisQuery(
            self,
            measurement_name,
            obs_query=obs_query or query.AxisQuery(),
            var_query=var_query or query.AxisQuery(),
            index_factory=functools.partial(
                IntIndexer,
                context=self.context,
            ),
        )

    def obs_axis_delete(
        self,
        coords: options.SparseDFCoords = (),
        *,
        value_filter: str | None = None,
        platform_config: options.PlatformConfig | None = None,
    ) -> None:
        """Delete observations (obs axis elements only) from all predefined sub-objects within an Experiment.

        This operation affects all Measurements within an Experiment. Deleted observations are specified:

        - by ``obs`` DataFrame soma_joinid: a sequence of soma_joinid or a slice
        - or, by ``obs`` DataFrame value filter / query condition

        Sub-objects affected by this operation:

        - exp.obs - delete matching DataFrame elements by soma_joinid
        - exp.obs_spatial_presence - delete matching DataFrame elements by soma_joinid
        - exp.ms[*].X[*] - delete matching SparseNDMatrix values, joining on the soma_dim_0 dimension
        - exp.ms[*].obsm[*] - delete matching SparseNDMatrix values, joining on the soma_dim_0 dimension
        - exp.ms[*].obsp[*] - delete matching SparseNDMatrix values, joining on the soma_dim_0 and soma_dim_1 dimension
        - exp.spatial[*].obsl[*] - delete matching DataFrame elements, joining on the soma_joinid dimension

        The following sub-objects are unaffected by this operation:

        - exp.ms[*].var
        - exp.ms[*].varm[*]
        - exp.ms[*].varp[*]
        - any user-defined arrays within the collection

        If any affected sub-objects are of type DenseNDArray (dense TileDB array), the entire operation results in an error.

        Either ``coords`` or ``value_filter`` must be provided. When both ``coords`` and ``value_filter`` are provided,
        the observations that match both constraints will be removed.

        For example, to delete observations where ``n_genes > 1000`` and ``n_counts < 2000``:
            >>> with tiledbsoma.Experiment.open(experiment_uri, mode="d") as exp:
            ...     exp.obs_axis_delete(value_filter="n_genes > 1000 and n_counts < 2000")

        Note: Deleting observations does not change the size of the current domain or possible enumeration values.

        Args:
            coords:
                A per-dimension ``Sequence`` of scalar, slice, sequence of scalar or
                `Arrow IntegerArray <https://arrow.apache.org/docs/python/generated/pyarrow.IntegerArray.html>` values
                defining the region to read.
            value_filter:
                An optional [value filter] to apply to the results.
                Defaults to no filter.
        """
        # build a list of arrays from which cells will be deleted.
        candidates = _create_obs_axis_candidates(self)

        # preemptive error checks
        if platform_config is not None and not isinstance(platform_config, TileDBDeleteOptions):
            raise TypeError(
                f"Invalid PlatformConfig with type {type(platform_config)}. Must have type {TileDBDeleteOptions.__name__}."
            )

        # query
        joinids = _query_joinids(self.obs.uri, coords, value_filter, self.context)

        # delete
        _delete_cells(candidates, joinids, platform_config=platform_config)

    def var_axis_delete(
        self,
        measurement_name: str,
        coords: options.SparseDFCoords = (),
        *,
        value_filter: str | None = None,
        platform_config: options.PlatformConfig | None = None,
    ) -> None:
        """Delete variables (``var`` axis elements only) from all predefined sub-objects within a single user-specified Measurement.

        Deleted features are specified:

        - by ``var`` DataFrame joinid: a sequence of joinid or a slice
        - by ``var`` DataFrame value filter / query condition

        Sub-objects affected by this operation, given a named Measurement (ms_name):

        - exp.ms[*ms_name*].var - delete matching DataFrame elements
        - exp.ms[*ms_name*].var_spatial_presence - delete matching DataFrame elements, joining on the soma_joinid dimension
        - exp.ms[*ms_name*].X[*] - delete matching SparseNDMatrix values, joining on the soma_dim_1 dimension
        - exp.ms[*ms_name*].varm[*] - delete matching SparseNDMatrix values, joining on the soma_dim_0 dimension
        - exp.ms[*ms_name*].varp[*] - delete matching SparseNDMatrix values, joining on the soma_dim_0 and soma_dim_1 dimensions
        - exp.spatial[*].varl[*ms_name*][*] - delete matching PointCloudDataFrame/GeometryDataFrame elements, joining on the soma_joinid dimension

        The following Measurement sub-objects are unaffected by this operation:

        - exp.ms[*ms_name*].obsm[*]
        - exp.ms[*ms_name*].obsp[*]

        If any affected sub-objects are a DenseNDArray (dense TileDB array), the entire operation results in an error.

        Either ``coords`` or ``value_filter`` must be provided. When both ``coords`` and ``value_filter`` are provided,
        the variables that match both constraints will be removed.

        For example, to delete observations where ``feature_biotype == 'gene'``:
            >>> with tiledbsoma.Experiment.open(experiment_uri, mode="d") as exp:
            ...     exp.var_axis_delete(value_filter="feature_biotype == 'gene'")

        Note: Deleting variables does not change the size of the current domain or possible enumeration values.

        Args:
            coords:
                A per-dimension ``Sequence`` of scalar, slice, sequence of scalar or
                `Arrow IntegerArray <https://arrow.apache.org/docs/python/generated/pyarrow.IntegerArray.html>` values
                defining the region to read.
            value_filter:
                An optional [value filter] to apply to the results.
                Defaults to no filter.
        """
        # build a list of arrays from which cells will be deleted.
        candidates = _create_var_axis_candidates(self, measurement_name)

        # preemptive error checks
        if platform_config is not None and not isinstance(platform_config, TileDBDeleteOptions):
            raise TypeError(
                f"Invalid PlatformConfig with type {type(platform_config)}. Must have type {TileDBDeleteOptions.__name__}."
            )

        # query
        joinids = _query_joinids(self.ms[measurement_name].var.uri, coords, value_filter, self.context)

        # delete
        _delete_cells(candidates, joinids, platform_config=platform_config)


@attrs.define(frozen=True)
class _ArrayDelMd:
    """Private helper class for axis deletes."""

    obj: DataFrame | SparseNDArray | PointCloudDataFrame | GeometryDataFrame
    name: str  # human-readable name
    join_on: tuple[str, ...]
    joinid_max: tuple[int, ...]


def _append_if_supported(candidates: list[_ArrayDelMd], obj: SOMAObject, name: str, join_on: tuple[str, ...]) -> None:
    """Private method.

    Append the candidate array to the final selection list if:

    - not a dense array
    - is one of the supported types
    - contains the expected join-on column.

    Primary goal here is to be graceful when the user has added ad hoc or otherwise
    surprising objects to the various Experiment sub-collections (e.g., obsm). While
    this should never occur, in practice it may.

    Dense array specifically causes an error. The other conditions cause the array
    in question to be ignored silently.
    """
    if isinstance(obj, DenseNDArray):
        raise SOMAError(f"Delete operation not permitted on dense array: {name}")
    if isinstance(obj, GeometryDataFrame):
        raise NotImplementedError(f"Delete operation not permitted on GeometryDataFrame: {name}")

    if not isinstance(obj, (DataFrame, SparseNDArray, PointCloudDataFrame)):
        return
    if not all(obj.schema.get_field_index(col_name) >= 0 for col_name in join_on):
        return

    if isinstance(obj, SparseNDArray):
        joinid_max = [obj.shape[int(col_name.removeprefix("soma_dim_"))] - 1 for col_name in join_on]
    else:
        # dataframe-esque
        joinid_max = [obj.domain[obj.index_column_names.index(col_name)][1] for col_name in join_on]

    candidates.append(_ArrayDelMd(obj=obj, name=name, join_on=join_on, joinid_max=tuple(joinid_max)))


def _create_obs_axis_candidates(exp: Experiment) -> list[_ArrayDelMd]:
    """Generate a list of candidate arrays for an obs axis delete."""
    arr: SOMAObject
    candidates: list[_ArrayDelMd] = []
    if "obs_spatial_presence" in exp:
        _append_if_supported(candidates, exp.obs_spatial_presence, "obs_spatial_presence", ("soma_joinid",))
    for ms_name, ms in exp.ms.items():
        if "X" in ms:
            for key, arr in ms.X.items():
                _append_if_supported(candidates, arr, f"ms[{ms_name}].X[{key}]", ("soma_dim_0",))
        if "obsm" in ms:
            for key, arr in ms.obsm.items():
                _append_if_supported(candidates, arr, f"ms[{ms_name}].obsm[{key}]", ("soma_dim_0",))
        if "obsp" in ms:
            for key, arr in ms.obsp.items():
                _append_if_supported(candidates, arr, f"ms[{ms_name}].obsp[{key}]", ("soma_dim_0", "soma_dim_1"))
    if "spatial" in exp:
        for scene_name, sc in exp.spatial.items():
            if "obsl" in sc:
                for key, arr in sc.obsl.items():
                    _append_if_supported(candidates, arr, f"spatial[{scene_name}].obsl[{key}]", ("soma_joinid",))

    _append_if_supported(candidates, exp.obs, "obs", ("soma_joinid",))
    return candidates


def _create_var_axis_candidates(exp: Experiment, ms_name: str) -> list[_ArrayDelMd]:
    """Generate a list of candidate arrays for a var axis delete."""
    arr: SOMAObject
    if ms_name not in exp.ms:
        raise ValueError(f"Measurement name {ms_name} does not exist in the experiment")

    candidates: list[_ArrayDelMd] = []
    if "var_spatial_presence" in exp.ms[ms_name]:
        _append_if_supported(
            candidates, exp.ms[ms_name].var_spatial_presence, f"ms[{ms_name}].var_spatial_presence", ("soma_joinid",)
        )
    if "X" in exp.ms[ms_name]:
        for key, arr in exp.ms[ms_name].X.items():
            _append_if_supported(candidates, arr, f"ms[{ms_name}].X[{key}]", ("soma_dim_1",))
    if "varm" in exp.ms[ms_name]:
        for key, arr in exp.ms[ms_name].varm.items():
            _append_if_supported(candidates, arr, f"ms[{ms_name}].varm[{key}]", ("soma_dim_0",))
    if "varp" in exp.ms[ms_name]:
        for key, arr in exp.ms[ms_name].varp.items():
            _append_if_supported(candidates, arr, f"ms[{ms_name}].varp[{key}]", ("soma_dim_0", "soma_dim_1"))
    if "spatial" in exp:
        for scene_name, sc in exp.spatial.items():
            if ("varl" in sc) and (ms_name in sc.varl):
                for key, arr in sc.varl[ms_name].items():
                    _append_if_supported(
                        candidates, arr, f"spatial[{scene_name}].varl[{ms_name}][{key}]", ("soma_joinid",)
                    )

    _append_if_supported(candidates, exp.ms[ms_name].var, f"ms[{ms_name}].var", ("soma_joinid",))
    return candidates


def _query_joinids(
    uri: str, coords: options.SparseDFCoords, value_filter: str | None, context: SOMAContext
) -> pa.Int64Array:
    with DataFrame.open(uri, context=context) as df:
        return (
            df
            .read(coords, value_filter=value_filter, column_names=["soma_joinid"])
            .concat()
            .column("soma_joinid")
            .combine_chunks()
        )


def _delete_cells(
    candidates: list[_ArrayDelMd], joinids: pa.Int64Array, platform_config: TileDBDeleteOptions | None = None
) -> None:
    if not len(joinids):
        return

    for arr_md in candidates:
        for join_on, joinid_max in zip(arr_md.join_on, arr_md.joinid_max):
            dim_idx = arr_md.obj.schema.get_field_index(join_on)
            coords = [slice(None)] * (dim_idx + 1)
            coords[dim_idx] = pa.compute.filter(joinids, pa.compute.less_equal(joinids, joinid_max))

            # Delete each join_on dimension separately, as we have no outer join
            arr_md.obj.delete_cells(tuple(coords), platform_config=platform_config)
