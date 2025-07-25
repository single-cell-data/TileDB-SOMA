# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

from __future__ import annotations

import warnings
from collections import defaultdict
from pathlib import Path
from typing import Any, Callable, Iterable, Sequence, cast

import anndata as ad
import attrs
import numpy as np
import numpy.typing as npt
import pandas as pd
import pyarrow as pa
from typing_extensions import Self, TypeAlias

import tiledbsoma
import tiledbsoma.io
import tiledbsoma.logging as logging
from tiledbsoma import DataFrame, Experiment, SOMAError
from tiledbsoma.options import SOMATileDBContext
from tiledbsoma.options._soma_tiledb_context import _validate_soma_tiledb_context

from .._util import read_h5ad
from .enum import extend_enumerations, get_enumerations
from .id_mappings import AxisIDMapping, ExperimentIDMapping, get_dataframe_values

MeasurementName: TypeAlias = str
ColumnName: TypeAlias = str


@attrs.define(kw_only=True, frozen=True)
class AxisAmbientLabelMapping:
    """For all the to-be-appended AnnData/H5AD inputs in SOMA multi-file append-mode ingestion, this
    class tracks the mapping of input-data ``obs`` or ``var`` ID-column name (barcode ID, gene
    symbol) to SOMA join IDs for SOMA experiment ``obs`` or ``var``, as well as any dictionary/enumeration
    values.

    See module-level comments for more information.
    """

    field_name: str
    joinid_map: pd.DataFrame  # id -> soma_joinid
    enum_values: dict[str, pd.CategoricalDtype]

    shape: int = attrs.field(init=False)

    def __attrs_post_init__(self) -> None:
        assert self.joinid_map.empty or self.joinid_map.soma_joinid.dtype == np.int64
        object.__setattr__(
            self,
            "shape",
            (
                int(self.joinid_map.soma_joinid.max() + 1)
                if len(self.joinid_map) > 0
                else 0
            ),
        )

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, AxisAmbientLabelMapping):
            raise NotImplementedError("Cannot compare to non-AxisAmbientLabelMapping")
        return (
            self.field_name == other.field_name
            and self.joinid_map.equals(other.joinid_map)
            and self.enum_values == other.enum_values
        )

    def id_mapping_from_values(self, input_ids: npt.ArrayLike) -> AxisIDMapping:
        """Given registered label-to-SOMA-join-ID mappings for all registered input files for an
        ``obs`` or ``var`` axis, and a list of input-file 0-up offsets, this returns an int-to-int
        mapping from a single input file's ``obs`` or ``var`` axis to the registered SOMA join IDs.
        """
        new_joinid_map = self.joinid_map.reindex(labels=input_ids, fill_value=-1)
        if new_joinid_map.soma_joinid.isin([-1]).any():
            raise ValueError(
                f"The input_ids for {self.field_name} [{input_ids[:10]}...] were not found in registration data."
            )
        return AxisIDMapping(data=new_joinid_map.soma_joinid.to_numpy())

    def id_mapping_from_dataframe(self, df: pd.DataFrame) -> AxisIDMapping:
        """Given registered label-to-SOMA-join-ID mappings for all registered input files for an
        ``obs`` or ``var`` axis, and an input file's dataframe with its 0-up offsets, this returns
        an int-to-int mapping from a single input file's ``obs`` or ``var`` axis to the registered
        SOMA join IDs.
        """
        values = get_dataframe_values(df, self.field_name)
        return self.id_mapping_from_values(values)


@attrs.define(kw_only=True, frozen=True)
class ExperimentAmbientLabelMapping:
    """For all the to-be-appended AnnData/H5AD inputs in SOMA multi-file append-mode ingestion, this
    class contains information required to perform ingestion via ``from_h5ad`` or ``from_anndata``.

    This class tracks the mapping from input-data ``obs`` or ``var`` ID-column name (barcode ID, gene
    symbol) to SOMA join IDs for SOMA experiment ``obs`` or ``var``, as well as any dictionary/enumeration
    values.
    """

    obs_axis: AxisAmbientLabelMapping
    var_axes: dict[str, AxisAmbientLabelMapping]
    prepared: bool = False

    def id_mappings_for_anndata(
        self, adata: ad.AnnData, *, measurement_name: str = "RNA"
    ) -> ExperimentIDMapping:

        obs_axis = AxisIDMapping(
            data=self.obs_axis.joinid_map.loc[
                _get_dataframe_joinid_index(adata.obs, self.obs_axis.field_name)
            ].soma_joinid.to_numpy()
        )
        var_axes = {
            measurement_name: AxisIDMapping(
                data=self.var_axes[measurement_name]
                .joinid_map.loc[
                    _get_dataframe_joinid_index(
                        adata.var, self.var_axes[measurement_name].field_name
                    )
                ]
                .soma_joinid.to_numpy()
            )
        }
        if adata.raw is not None:
            var_axes["raw"] = AxisIDMapping(
                data=self.var_axes["raw"]
                .joinid_map.loc[
                    _get_dataframe_joinid_index(
                        adata.raw.var, self.var_axes["raw"].field_name
                    )
                ]
                .soma_joinid.to_numpy()
            )

        return ExperimentIDMapping(obs_axis=obs_axis, var_axes=var_axes)

    def get_obs_shape(self) -> int:
        return self.obs_axis.shape

    def get_var_shapes(self) -> dict[str, int]:
        return {ms_name: self.var_axes[ms_name].shape for ms_name in self.var_axes}

    def subset_for_anndata(self, adata: ad.AnnData) -> Self:
        """Return a copy of this object containing only the information necessary to ingest
        the specified AnnData.

        This is an optional step, used to improve the performance of multi-process or
        distributed ingestion. It reduces the registration information transmitted to
        the worker process performing the call to ``from_anndata`` or ``from_h5ad``.

        Args:
            h5ad_path: an ``anndata.AnnData``, previously registered in this ``ExperimentAmbientLabelMapping``.

        Returns:
            A new ``ExperimentAmbientLabelMapping`` scoped specifically for the AnnData.
        """
        # Only subset obs - provides largest benefit with simple implementation.
        return attrs.evolve(
            self,
            obs_axis=attrs.evolve(
                self.obs_axis,
                joinid_map=self.obs_axis.joinid_map.loc[
                    _get_dataframe_joinid_index(adata.obs, self.obs_axis.field_name)
                ],
            ),
        )

    def subset_for_h5ad(self, h5ad_path: str) -> Self:
        """Subset this plan to only contain ID maps useful for this H5AD.

        See ``subset_for_anndata`` for more information.

        Args:
            h5ad_path: path to H5AD

        Returns:
            A new ``ExperimentAmbientLabelMapping`` scoped specifically for the H5AD.
        """
        with read_h5ad(h5ad_path, mode="r") as adata:
            return self.subset_for_anndata(adata)

    def prepare_experiment(
        self, experiment_uri: str, context: SOMATileDBContext | None = None
    ) -> None:
        """Prepare experiment for ingestion.

        Currently performs two operations:
        1. Resize experiment to a shape sufficient to contain all registered AnnData/H5AD inputs
        2. Evolve schema on all dict/enum/categorical columns to include any new values defined in registered AnnData (e.g., Pandas Categoricals with additional categories).

        This makes subsequent data writes race safe, for workflows using concurrent dataset writers
        (ie., parallel calls to `to_anndata` or `from_h5ad`).

        This operation must be performed after the experiment is created, and before any writes
        to the experiment.

        Args:
            experiment_uri: the Experiment to prepare for ingestion.

            context: a SOMA context

        Returns:
            None
        """
        context = _validate_soma_tiledb_context(context)

        def _check_experiment_structure(exp: tiledbsoma.Experiment) -> None:
            # Verify that the experiment has been created correctly - check for existence of obs & var
            # and raise error if Experiment does not contain expected structure.
            did_you_create = (
                "Did you create the Experiment using `from_anndata` or `from_h5ad`?"
            )
            if "obs" not in exp:
                raise ValueError(
                    f"SOMA Experiment is missing required 'obs' DataFrame. {did_you_create}"
                )
            for ms_name in self.var_axes:
                if len(self.var_axes[ms_name].joinid_map) > 0:
                    if ms_name not in exp.ms:
                        raise ValueError(
                            f"SOMA Experiment is missing required Measurement '{ms_name}'. {did_you_create}"
                        )
                    if "var" not in exp.ms[ms_name]:
                        raise ValueError(
                            f"SOMA Experiment is missing required `var` Dataframe in Measurement '{ms_name}'. {did_you_create}"
                        )

        with Experiment.open(experiment_uri, context=context) as E:

            _check_experiment_structure(E)

            # Resize is done only if we have an Experiment supporting current domain.
            # Code assumes that if obs is of a given era (i.e. pre/post current domain change),
            # so are all other arrays in the experiment.
            ok_to_resize, _ = E.obs.tiledbsoma_resize_soma_joinid_shape(
                self.get_obs_shape(), check_only=True
            )
            if ok_to_resize:
                tiledbsoma.io.resize_experiment(
                    experiment_uri,
                    nobs=self.get_obs_shape(),
                    nvars=self.get_var_shapes(),
                    context=context,
                )
            else:
                warnings.warn(
                    "Experiment does not support resizing. Please consider upgrading the dataset "
                    "using 'tiledbsoma.io.upgrade_experiment_shapes'."
                )

        with Experiment.open(experiment_uri, context=context, mode="w") as E:

            # Enumerations schema evolution
            extend_enumerations(E.obs, self.obs_axis.enum_values)

            for ms_name, var_axis in self.var_axes.items():
                if var_axis.enum_values:
                    extend_enumerations(E.ms[ms_name].var, var_axis.enum_values)

        # The class is a frozen `attrs` instance, to protect from user modification of the data.
        # This is the "blessed" way for an implementation to modify itself (per attrs docs).
        object.__setattr__(self, "prepared", True)

    @staticmethod
    def _load_axes_metadata_from_anndatas(
        adatas: Iterable[ad.AnnData],
        obs_field_name: str,
        var_field_name: str,
        validate_anndata: Callable[[ad.AnnData], None],
    ) -> tuple[AnnDataAxisMetadata, AnnDataAxisMetadata, AnnDataAxisMetadata | None]:
        """Private helper to load axis metadata for obs, var and raw.var (if present) from AnnData,
        in an intermediate form.

        Return (obs, var, raw.var)
        """
        obs_metadata: list[AnnDataAxisMetadata] = []
        var_metadata: list[AnnDataAxisMetadata] = []
        raw_var_metadata: list[AnnDataAxisMetadata] = []

        def categorical_columns(df: pd.DataFrame) -> dict[Any, pd.CategoricalDtype]:
            return cast(
                dict[str, pd.CategoricalDtype],
                {k: v.dtype for k, v in df.items() if v.dtype == "category"},
            )

        for adata in adatas:

            validate_anndata(adata)  # may throw

            obs_metadata.append(
                AnnDataAxisMetadata(
                    field_name=obs_field_name,
                    field_index=_get_dataframe_joinid_index(adata.obs, obs_field_name),
                    enum_values=categorical_columns(adata.obs),
                )
            )
            var_metadata.append(
                AnnDataAxisMetadata(
                    field_name=var_field_name,
                    field_index=_get_dataframe_joinid_index(adata.var, var_field_name),
                    enum_values=categorical_columns(adata.var),
                )
            )
            if adata.raw is not None:
                raw_var_metadata.append(
                    AnnDataAxisMetadata(
                        field_name=var_field_name,
                        field_index=_get_dataframe_joinid_index(
                            adata.raw.var, var_field_name
                        ),
                        enum_values=categorical_columns(adata.raw.var),
                    )
                )

        obs, var, raw_var = (
            AnnDataAxisMetadata.reduce(obs_metadata),
            AnnDataAxisMetadata.reduce(var_metadata),
            AnnDataAxisMetadata.reduce(raw_var_metadata),
        )
        assert obs is not None
        assert var is not None
        return obs, var, raw_var

    @staticmethod
    def _load_axes_metadata_from_h5ads(
        paths: Sequence[str | Path],
        obs_field_name: str,
        var_field_name: str,
        validate_anndata: Callable[[ad.AnnData], None],
    ) -> tuple[AnnDataAxisMetadata, AnnDataAxisMetadata, AnnDataAxisMetadata | None]:
        """Private helper to load axis metadata for obs, var and raw.var (if present) from H5ADs.

        Return (obs, var, raw.var).
        """
        obs_metadata: list[AnnDataAxisMetadata] = []
        var_metadata: list[AnnDataAxisMetadata] = []
        raw_var_metadata: list[AnnDataAxisMetadata] = []

        for p in paths:
            with read_h5ad(p, mode="r") as adata:
                obs, var, raw_var = (
                    ExperimentAmbientLabelMapping._load_axes_metadata_from_anndatas(
                        [adata], obs_field_name, var_field_name, validate_anndata
                    )
                )
            obs_metadata.append(obs)
            var_metadata.append(var)
            if raw_var is not None:
                raw_var_metadata.append(raw_var)

        obs, var, raw_var = (
            AnnDataAxisMetadata.reduce(obs_metadata),
            AnnDataAxisMetadata.reduce(var_metadata),
            AnnDataAxisMetadata.reduce(raw_var_metadata),
        )
        assert obs is not None
        assert var is not None
        return obs, var, raw_var

    @staticmethod
    def _load_existing_experiment_metadata(
        uri: str, obs_field_name: str, var_field_name: str, context: SOMATileDBContext
    ) -> tuple[
        pd.DataFrame,
        dict[ColumnName, pd.DataFrame],
        dict[ColumnName, pd.CategoricalDtype],
        dict[MeasurementName, dict[ColumnName, pd.CategoricalDtype]],
    ]:
        """Private helper to load any joinid/enum metadata from an existing experiment.

        Returns (obs_joinid_map, var_joinid_maps, obs_enum_values, var_enum_values).
        """

        def _get_enum_values(df: DataFrame) -> dict[str, pd.CategoricalDtype]:
            return get_enumerations(
                df, [f.name for f in df.schema if pa.types.is_dictionary(f.type)]
            )

        def _get_joinid_map(df: DataFrame, field_name: str) -> pd.DataFrame:
            return cast(
                pd.DataFrame,
                df.read(column_names=["soma_joinid", field_name])
                .concat()
                .to_pandas()
                .set_index(field_name),
            )

        with Experiment.open(uri, context=context) as E:
            existing_obs_joinid_map = _get_joinid_map(E.obs, obs_field_name)
            existing_obs_enum_values = _get_enum_values(E.obs)
            existing_var_joinid_maps = {}
            existing_var_enum_values = {}
            for ms_name in E.ms.keys():
                if (
                    "var" in E.ms[ms_name]
                    and var_field_name in E.ms[ms_name].var.keys()
                ):
                    existing_var_joinid_maps[ms_name] = _get_joinid_map(
                        E.ms[ms_name].var, var_field_name
                    )
                    existing_var_enum_values[ms_name] = _get_enum_values(
                        E.ms[ms_name].var
                    )

        return (
            existing_obs_joinid_map,
            existing_var_joinid_maps,
            existing_obs_enum_values,
            existing_var_enum_values,
        )

    @staticmethod
    def _register_common(
        experiment_uri: str | None,
        axes_metadata: list[
            tuple[AnnDataAxisMetadata, AnnDataAxisMetadata, AnnDataAxisMetadata | None]
        ],
        *,
        measurement_name: str,
        obs_field_name: str,
        var_field_name: str,
        context: SOMATileDBContext,
    ) -> ExperimentAmbientLabelMapping:
        """Private method used by various constructor paths -- shared code for registration.

        Four-step process common to all ingestion plans. Given AnnData axis metadata (joinid map
        and enum definitions):
        1. Load plan info from existing Experiment, if available
        2. Reduce all AnnData axis metadata
        3. Create a joinid map for ingestion
        4. Create enum evolution
        5. Create and return the plan.

        [sc-65318]: the current design assumes that the join column for `var` is the same across
        all measurements. For simple cases (e.g., single modality) this is normally true.
        It is less likely to be true for multi-modal data. Future rework should address this
        by allowing a per-measurement var join column.
        """
        tp = context.threadpool
        existing_obs_joinid_map: pd.DataFrame
        existing_var_joinid_maps: dict[str, pd.DataFrame]
        existing_obs_enum_values: dict[str, pd.CategoricalDtype]
        existing_var_enum_values: dict[str, dict[str, pd.CategoricalDtype]]

        #
        # Step 1: load all existing Experiment info.
        #
        if experiment_uri is not None:
            logging.log_io_same("Loading existing experiment joinid map")
            experiment_metadata_ft = tp.submit(
                ExperimentAmbientLabelMapping._load_existing_experiment_metadata,
                experiment_uri,
                obs_field_name,
                var_field_name,
                context,
            )

        #
        # Step 2: reduce axis metadata
        #
        logging.log_io_same("Reducing axis metadata")
        obs_axis_metadata, var_axis_metadata, raw_var_axis_metadata = tp.map(
            AnnDataAxisMetadata.reduce,
            [
                [t[0] for t in axes_metadata if t[0] is not None],  # obs
                [t[1] for t in axes_metadata if t[1] is not None],  # var
                [t[2] for t in axes_metadata if t[2] is not None],  # raw.var
            ],
        )
        logging.log_io_same("Finished reducing axis metadata")

        # And, grab the result of step 1 from the futures
        if experiment_uri is not None:
            (
                existing_obs_joinid_map,
                existing_var_joinid_maps,
                existing_obs_enum_values,
                existing_var_enum_values,
            ) = experiment_metadata_ft.result()
            logging.log_io_same("Existing joinid maps are loaded.")
        else:
            existing_obs_joinid_map = pd.DataFrame()
            existing_var_joinid_maps = {
                measurement_name: pd.DataFrame(),
                "raw": pd.DataFrame(),
            }
            existing_obs_enum_values = {}
            existing_var_enum_values = {measurement_name: {}, "raw": {}}

        #
        # Step 3: create a joinid map for each axis
        #
        def _make_joinid_map(
            joinids_index: pd.Index,
            prev_joinid_map: pd.DataFrame,
        ) -> pd.DataFrame:
            maps = []
            if len(prev_joinid_map) > 0:
                joinids_index = joinids_index.difference(
                    prev_joinid_map.index, sort=False
                )
                next_soma_joinid = prev_joinid_map.soma_joinid.max() + 1
                maps.append(prev_joinid_map)
            else:
                next_soma_joinid = 0

            logging.log_io_same(f"next soma_joinid={next_soma_joinid}")
            maps.append(
                pd.DataFrame(
                    index=joinids_index,
                    data={
                        "soma_joinid": np.arange(
                            next_soma_joinid,
                            next_soma_joinid + len(joinids_index),
                            dtype=np.int64,
                        )
                    },
                )
            )
            return pd.concat(maps)

        obs_joinid_map_future = tp.submit(
            _make_joinid_map, obs_axis_metadata.field_index, existing_obs_joinid_map
        )
        var_joinid_maps_future = {
            measurement_name: tp.submit(
                _make_joinid_map,
                var_axis_metadata.field_index,
                existing_var_joinid_maps.get(measurement_name, pd.DataFrame()),
            )
        }
        if len(raw_var_axis_metadata.field_index) > 0:
            var_joinid_maps_future["raw"] = tp.submit(
                _make_joinid_map,
                raw_var_axis_metadata.field_index,
                existing_var_joinid_maps.get("raw", pd.DataFrame()),
            )

        obs_joinid_map = obs_joinid_map_future.result()
        var_joinid_maps = existing_var_joinid_maps | {
            k: f.result() for k, f in var_joinid_maps_future.items()
        }

        #
        # Step 4: create merged enum values for all axis dataframes
        #
        obs_enum_values = AnnDataAxisMetadata.reduce_enum_values(
            [existing_obs_enum_values, obs_axis_metadata.enum_values]
        )
        var_enum_values = existing_var_enum_values.copy()
        var_enum_values[measurement_name] = AnnDataAxisMetadata.reduce_enum_values(
            [
                existing_var_enum_values.get(measurement_name, {}),
                var_axis_metadata.enum_values,
            ]
        )
        if len(raw_var_axis_metadata.enum_values) > 0:
            var_enum_values["raw"] = AnnDataAxisMetadata.reduce_enum_values(
                [
                    existing_var_enum_values.get("raw", {}),
                    raw_var_axis_metadata.enum_values,
                ]
            )

        #
        # Step 5: create and return the ingestion plan
        #
        obs_axis = AxisAmbientLabelMapping(
            field_name=obs_field_name,
            joinid_map=obs_joinid_map,
            enum_values=obs_enum_values,
        )
        var_axes = {
            k: AxisAmbientLabelMapping(
                field_name=var_field_name,
                joinid_map=(
                    var_joinid_maps[k] if k in var_joinid_maps else pd.DataFrame()
                ),
                enum_values=var_enum_values[k] if k in var_enum_values else {},
            )
            for k in set(var_joinid_maps.keys()) | set(var_enum_values.keys())
        }
        return ExperimentAmbientLabelMapping(obs_axis=obs_axis, var_axes=var_axes)

    @staticmethod
    def _validate_anndata(append_obsm_varm: bool, adata: ad.AnnData) -> None:
        """Pre-checks performed on all AnnData."""

        def check_df(df: pd.DataFrame | None, df_name: str) -> None:
            if df is None or df.index.empty:
                raise ValueError(
                    f"Unable to ingest AnnData with empty {df_name} dataframe."
                )
            elif not df.index.is_unique:
                raise ValueError(
                    f"Non-unique registration values have been provided in {df_name} dataframe."
                )

        check_df(adata.obs, "obs")
        check_df(adata.var, "var")
        if adata.raw is not None:
            check_df(adata.raw.var, "raw.var")

        if not append_obsm_varm:
            if len(adata.obsm) > 0 or len(adata.varm) > 0:
                raise ValueError(
                    "The append-mode ingest of obsm and varm is only supported via explicit opt-in. Please drop them from the inputs, or retry with append_obsm_varm=True."
                )

        if len(adata.obsp) > 0 or len(adata.varp) > 0:
            raise ValueError(
                "The append-mode ingest of obsp and varp is not supported. Please retry without them."
            )

        if adata.uns:
            warnings.warn(
                "The append-mode ingest of 'uns' is typically an error due to uns key collisions "
                "across multiple AnnData. Drop 'uns' from AnnData to remove this warning, or if you "
                "intend for 'uns' to merge, ensure each AnnData uses unique keys."
            )


@attrs.define(kw_only=True, frozen=True)
class AnnDataAxisMetadata:
    """Private class encapsulating _intermediate_ information extracted from registered
    H5ADs/AnnData.

    Other than data storage, this also provides a reducer for the type.
    """

    field_name: str  # user-specified join field name
    field_index: pd.Index[Any]  # index of join field values
    enum_values: dict[ColumnName, pd.CategoricalDtype]

    @classmethod
    def reduce(cls, ams: list[Self]) -> AnnDataAxisMetadata:
        assert all(isinstance(a, cls) for a in ams)
        assert all([a.field_name == ams[0].field_name for a in ams])
        if not all(a.field_index.dtype == ams[0].field_index.dtype for a in ams):
            raise SOMAError("All AnnData must have a common dtype for their index.")

        if len(ams) == 0:
            return cls(field_name="", field_index=pd.Index([]), enum_values={})
        if len(ams) == 1:
            return ams[0]

        def _reduce_field_index(indices: list[pd.Index]) -> pd.Index:
            """Reducer for joinid indices."""
            if len(indices) == 0:
                return pd.Index([])
            if len(indices) == 1:
                return indices[0]
            return cast("pd.Index[Any]", indices[0].append(indices[1:]).drop_duplicates())  # type: ignore[no-untyped-call]

        return cls(
            field_name=ams[0].field_name,
            field_index=_reduce_field_index(
                [a.field_index for a in ams if not a.field_index.empty]
            ),
            enum_values=cls.reduce_enum_values([a.enum_values for a in ams]),
        )

    @staticmethod
    def reduce_enum_values(
        enum_values: list[dict[ColumnName, pd.CategoricalDtype]],
    ) -> dict[ColumnName, pd.CategoricalDtype]:
        """Reducer for enum values."""

        def _merge_categoricals(
            col_name: str,
            enums: list[pd.CategoricalDtype],
        ) -> pd.CategoricalDtype:
            assert len(enums) > 0
            if len(enums) == 1:
                return enums[0]

            ordered = enums[0].ordered
            if not all(e.ordered == ordered for e in enums[1:]):
                raise SOMAError(
                    f"Unable to register AnnData -- for column `{col_name}`, all AnnData dtype must have the same categorical ordering."
                )

            if not ordered:
                return pd.CategoricalDtype(
                    enums[0]
                    .categories.append([e.categories for e in enums[1:]])  # type: ignore[no-untyped-call]
                    .drop_duplicates(),
                    ordered=False,
                )

            # Ordered enums are tricky - handle the simple case where all are identical
            # and error on anything else.
            for e in enums[1:]:
                if e != enums[0]:
                    raise SOMAError(
                        f"Unable to register AnnData -- for column `{col_name}`, all AnnData must have the same dtype."
                    )
            return enums[0]

        # invert the enum maps from list[dict[str, pd.Categorical]] to dict[str, list[pd.Categorical]]
        inverted_enum_values: dict[str, list[pd.CategoricalDtype]] = defaultdict(list)
        for e in enum_values:
            for k, v in e.items():
                inverted_enum_values[k].append(v)

        return {k: _merge_categoricals(k, v) for k, v in inverted_enum_values.items()}


def _get_dataframe_joinid_index(df: pd.DataFrame, field_name: str) -> pd.Index:
    """Given an AnnData obs/var, extract the index for the user-selected join column."""
    if field_name in df:
        return cast("pd.Index[Any]", pd.Index(df[field_name]))
    if df.index.name in (field_name, "index", None):
        return cast("pd.Index[Any]", df.index)
    raise ValueError(f"Could not find field name {field_name} in dataframe.")
