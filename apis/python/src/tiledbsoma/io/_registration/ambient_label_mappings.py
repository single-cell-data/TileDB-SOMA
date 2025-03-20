# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

from __future__ import annotations

import json
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
from typing_extensions import Self

import tiledbsoma
import tiledbsoma.logging as logging
from tiledbsoma import DataFrame, Experiment, SOMAError
from tiledbsoma.options import SOMATileDBContext

from .enum import _extend_enumeration, _get_enumeration
from .id_mappings import AxisIDMapping, ExperimentIDMapping, get_dataframe_values


@attrs.define(kw_only=True, frozen=True)
class AxisAmbientLabelMapping:
    """
    For all the to-be-appended AnnData/H5AD inputs in SOMA multi-file append-mode ingestion, this
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

    def get_shape(self) -> int:
        return self.shape

    def id_mapping_from_values(self, input_ids: npt.ArrayLike) -> AxisIDMapping:
        """Given registered label-to-SOMA-join-ID mappings for all registered input files for an
        ``obs`` or ``var`` axis, and a list of input-file 0-up offsets, this returns an int-to-int
        mapping from a single input file's ``obs`` or ``var`` axis to the registered SOMA join IDs.
        """
        new_joinid_map = self.joinid_map.reindex(labels=input_ids, fill_value=-1)
        if new_joinid_map.soma_joinid.isin([-1]).any():
            raise ValueError(
                f"input_ids for {self.field_name} [{input_ids[:10]}...] not found in registration data"
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

    @staticmethod
    def _attrs_todict_serializer(
        _: Any, field: attrs.Attribute | None, value: Any  # type: ignore[type-arg]
    ) -> Any:
        if field is None:
            return value
        if field.name == "joinid_map":
            return value.to_dict(orient="tight")
        if field.name == "enum_values":
            return {
                k: dict(ordered=v.ordered, categories=v.categories.tolist())
                for k, v in value.items()
            }
        if field.name == "shape":
            return int(value)
        if field.name == "field_name":
            return value

        # else, unknown field, raise
        raise NotImplementedError("Unknown field - should never happen.")

    @staticmethod
    def _attrs_fromdict_deserializer(field: str, value: Any) -> Any:
        if field == "joinid_map":
            return pd.DataFrame.from_dict(value, orient="tight")
        if field == "enum_values":
            return {
                k: pd.CategoricalDtype(ordered=v["ordered"], categories=v["categories"])
                for k, v in value.items()
            }
        return value

    @classmethod
    def from_dict(cls, d: dict[Any, Any]) -> Self:
        return cls(
            **{
                k: cls._attrs_fromdict_deserializer(k, v)
                for k, v in d.items()
                if k not in ["shape"]  # handled in post-init
            }
        )

    def to_json(self) -> str:
        return json.dumps(
            attrs.asdict(self, value_serializer=self._attrs_todict_serializer)
        )

    @classmethod
    def from_json(cls, s: str) -> Self:
        return cls.from_dict(json.loads(s))


@attrs.define(kw_only=True, frozen=True)
class ExperimentAmbientLabelMapping:
    """
    For all the to-be-appended AnnData/H5AD inputs in SOMA multi-file append-mode ingestion, this
    class contains information required to perform ingestion via ``from_h5ad`` or ``from_anndata``.

    This class tracks the mapping from input-data ``obs`` or ``var`` ID-column name (barcode ID, gene
    symbol) to SOMA join IDs for SOMA experiment ``obs`` or ``var``, as well as any dictionary/enumeration
    values.
    """

    obs_axis: AxisAmbientLabelMapping
    var_axes: dict[str, AxisAmbientLabelMapping]

    def id_mappings_for_anndata(
        self, adata: ad.AnnData, *, measurement_name: str = "RNA"
    ) -> ExperimentIDMapping:

        obs_axis = AxisIDMapping(
            data=self.obs_axis.joinid_map.loc[adata.obs.index].soma_joinid.to_numpy()
        )
        var_axes = {
            measurement_name: AxisIDMapping(
                data=self.var_axes[measurement_name]
                .joinid_map.loc[adata.var.index]
                .soma_joinid.to_numpy()
            )
        }
        if adata.raw is not None:
            var_axes["raw"] = AxisIDMapping(
                data=self.var_axes["raw"]
                .joinid_map.loc[adata.raw.var.index]
                .soma_joinid.to_numpy()
            )

        return ExperimentIDMapping(obs_axis=obs_axis, var_axes=var_axes)

    def get_obs_shape(self) -> int:
        return self.obs_axis.shape

    def get_var_shapes(self) -> dict[str, int]:
        return {ms_name: self.var_axes[ms_name].shape for ms_name in self.var_axes}

    def get_obs_enum_values(self) -> dict[str, pd.CategoricalDtype]:
        return self.obs_axis.enum_values

    def get_var_enum_values(self) -> dict[str, dict[str, pd.CategoricalDtype]]:
        return {
            ms_name: self.var_axes[ms_name].enum_values for ms_name in self.var_axes
        }

    def subset_for_anndata(self, adata: ad.AnnData) -> Self:
        """Return a copy of this object containing only the information necessary to ingest
        the specified AnnData.

        This is an optional step, used to improve the performance of multi-process or
        distributed ingestion. It reduces the registration information transmitted to
        the worker process performing the call to ``from_anndata`` or ``from_h5ad``.
        """

        # Just do obs - provides largest benefit with simple implementation.
        return attrs.evolve(
            self,
            obs_axis=attrs.evolve(
                self.obs_axis, joinid_map=self.obs_axis.joinid_map.loc[adata.obs.index]
            ),
        )

    def subset_for_h5ad(self, h5ad_path: str) -> Self:
        """Subset this plan to only contain ID maps useful for this H5AD. See ``subset_for_anndata``
        for more information."""
        adata = ad.read_h5ad(h5ad_path, backed="r")
        return self.subset_for_anndata(adata)

    def prepare_experiment(
        self, experiment_uri: str, context: SOMATileDBContext
    ) -> None:
        """Prepare experiment for ingestion.

        Currently performs two operations:
        1. Resize experiment to a shape sufficient to contain all registered AnnData
        2. Evolve schema on all dict/enum/categorical columns to include any new values defined in
           registered AnnData (e.g., Pandas Categoricals with additional categories).

        This operation must be performed before any writes to the experiment.
        """
        tiledbsoma.io.resize_experiment(
            experiment_uri, nobs=self.get_obs_shape(), nvars=self.get_var_shapes()
        )

        with Experiment.open(experiment_uri, context=context) as E:
            for k, v in self.obs_axis.enum_values.items():
                _extend_enumeration(E.obs, k, pa.array(v.categories))

            for ms_name, var_axis in self.var_axes.items():
                for k, v in var_axis.enum_values.items():
                    _extend_enumeration(E.ms[ms_name].var, k, pa.array(v.categories))

    def to_json(self) -> str:
        """The ``to_json`` and ``from_json`` methods allow you to persist
        the registration mappings to disk.

        Hint: pickle may be more efficient for some use cases."""
        dikt = {
            "obs_axis": attrs.asdict(
                self.obs_axis,
                recurse=True,
                value_serializer=AxisAmbientLabelMapping._attrs_todict_serializer,
            ),
            "var_axes": {
                k: attrs.asdict(
                    v,
                    recurse=True,
                    value_serializer=AxisAmbientLabelMapping._attrs_todict_serializer,
                )
                for k, v in self.var_axes.items()
            },
        }
        return json.dumps(dikt)

    @classmethod
    def from_json(cls, s: str) -> Self:
        """The ``to_json`` and ``from_json`` methods allow you to persist
        the registration mappings to disk."""
        d = json.loads(s)
        return cls(
            obs_axis=AxisAmbientLabelMapping.from_dict(d["obs_axis"]),
            var_axes={
                k: AxisAmbientLabelMapping.from_dict(v)
                for k, v in d["var_axes"].items()
            },
        )

    @staticmethod
    def _load_axes_metadata_from_anndatas(
        adatas: Iterable[ad.AnnData],
        obs_field_name: str,
        var_field_name: str,
        validate_anndata: Callable[[ad.AnnData], None],
    ) -> tuple[AxisMetadata, AxisMetadata, AxisMetadata | None]:
        """Load axis metadata for obs, var and raw.var (if present) from AnnData,
        in an intermediate form.

        Return (obs, var, raw.var)
        """
        obs_metadata: list[AxisMetadata] = []
        var_metadata: list[AxisMetadata] = []
        raw_var_metadata: list[AxisMetadata] = []

        for adata in adatas:

            validate_anndata(adata)  # may throw

            obs_metadata.append(
                AxisMetadata(
                    field_name=obs_field_name,
                    field_index=adata.obs.index,
                    enum_values={
                        k: v.dtype
                        for k, v in adata.obs.items()
                        if v.dtype == "category"
                    },
                )
            )
            var_metadata.append(
                AxisMetadata(
                    field_name=var_field_name,
                    field_index=adata.var.index,
                    enum_values={
                        k: v.dtype
                        for k, v in adata.var.items()
                        if v.dtype == "category"
                    },
                )
            )
            if adata.raw is not None:
                raw_var_metadata.append(
                    AxisMetadata(
                        field_name=var_field_name,
                        field_index=adata.raw.var.index,
                        enum_values={
                            k: v.dtype
                            for k, v in adata.raw.var.items()
                            if v.dtype == "category"
                        },
                    )
                )

        obs, var, raw_var = (
            AxisMetadata.reduce(obs_metadata),
            AxisMetadata.reduce(var_metadata),
            AxisMetadata.reduce(raw_var_metadata),
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
    ) -> tuple[AxisMetadata, AxisMetadata, AxisMetadata | None]:
        """Load axis metadata for obs, var and raw.var (if present) from H5ADs.

        Return (obs, var, raw.var)
        """
        return ExperimentAmbientLabelMapping._load_axes_metadata_from_anndatas(
            (ad.read_h5ad(path, backed="r") for path in paths),  # lazy open
            obs_field_name,
            var_field_name,
            validate_anndata,
        )

    @staticmethod
    def _load_existing_experiment_metadata(
        uri: str, obs_field_name: str, var_field_name: str, context: SOMATileDBContext
    ) -> tuple[
        pd.DataFrame,
        dict[str, pd.DataFrame],
        dict[str, pd.CategoricalDtype],  # { col_name: CategoricalDtype, ...}
        dict[
            str, dict[str, pd.CategoricalDtype]
        ],  # { ms_name: { col_name: CategoricalDtype, ...}, ...}
    ]:
        """Private helper to load any joinid/enum metadata from an existing experiment."""

        def _get_enum_values(df: DataFrame) -> dict[str, pd.CategoricalDtype]:
            return {
                f.name: _get_enumeration(df, f.name)
                for f in df.schema
                if pa.types.is_dictionary(f.type)
            }

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
        axes_metadata: list[tuple[AxisMetadata, AxisMetadata, AxisMetadata | None]],
        *,
        measurement_name: str,
        obs_field_name: str,
        var_field_name: str,
        context: SOMATileDBContext,
    ) -> ExperimentAmbientLabelMapping:
        """
        Private method used by various constructor paths.

        Four-step process common to all ingestion plans. Given AnnData axis metadata (joinid map
        and enum definitions):
        1. Load plan info from existing Experiment, if available
        2. Reduce all AnnData axis metadata
        3. Create a joinid map for ingestion
        4. Create enum evolution
        and return the plan.
        """

        tp = context.threadpool
        existing_obs_joinid_map: pd.DataFrame
        existing_var_joinid_maps: dict[str, pd.DataFrame]
        existing_obs_enum_values: dict[str, pd.CategoricalDtype]
        existing_var_enum_values: dict[str, dict[str, pd.CategoricalDtype]]

        #
        # Step 1: load all existing Experiment info
        #
        if experiment_uri is not None:
            logging.log_io(None, "Loading existing experiment joinid map")
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
        logging.log_io(None, "Reducing axis metadata")
        obs_axis_metadata, var_axis_metadata, raw_var_axis_metadata = tp.map(
            AxisMetadata.reduce,
            [
                [t[0] for t in axes_metadata if t[0] is not None],  # obs
                [t[1] for t in axes_metadata if t[1] is not None],  # var
                [t[2] for t in axes_metadata if t[2] is not None],  # raw.var
            ],
        )
        logging.log_io(None, "Finished reducing axis metadata")

        # And, grab the result of step 1
        if experiment_uri is not None:
            (
                existing_obs_joinid_map,
                existing_var_joinid_maps,
                existing_obs_enum_values,
                existing_var_enum_values,
            ) = experiment_metadata_ft.result()
            logging.log_io(None, "Existing joinid maps are loaded.")
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
            joinids_index: pd.Index,  # type:ignore[type-arg]
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

            logging.log_io(None, f"next soma_joinid={next_soma_joinid}")
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

        obs_enum_values = AxisMetadata.reduce_enum_values(
            [existing_obs_enum_values, obs_axis_metadata.enum_values]
        )
        var_enum_values = existing_var_enum_values.copy()
        var_enum_values[measurement_name] = AxisMetadata.reduce_enum_values(
            [
                existing_var_enum_values.get(measurement_name, {}),
                var_axis_metadata.enum_values,
            ]
        )
        if len(raw_var_axis_metadata.enum_values) > 0:
            var_enum_values["raw"] = AxisMetadata.reduce_enum_values(
                [
                    existing_var_enum_values.get("raw", {}),
                    raw_var_axis_metadata.enum_values,
                ]
            )

        # and return the ingestion plan
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
    def _validate_anndata(
        append_obsm_varm: bool,
        adata: ad.AnnData,
    ) -> None:
        """Pre-checks performed on all AnnData"""

        def check_df(df: pd.DataFrame | None, df_name: str) -> None:
            if df is None or df.empty:
                raise ValueError(
                    f"Unable to ingest AnnData with empty {df_name} dataframe"
                )
            elif not df.index.is_unique:
                raise ValueError(
                    f"non-unique registration values have been provided in {df_name} dataframe"
                )

        check_df(adata.obs, "obs")
        check_df(adata.var, "var")
        if adata.raw is not None:
            check_df(adata.raw.var, "raw.var")

        if not append_obsm_varm:
            if len(adata.obsm) > 0 or len(adata.varm) > 0:
                raise ValueError(
                    "append-mode ingest of obsm and varm is only supported via explicit opt-in. Please drop them from the inputs, or retry with append_obsm_varm=True."
                )

        if len(adata.obsp) > 0 or len(adata.varp) > 0:
            raise ValueError(
                "append-mode ingest of obsp and varp is not supported. Please retry without them."
            )

        if adata.uns:
            warnings.warn(
                "append-mode ingest of 'uns' is typically an error due to uns key collisions "
                "across multiple AnnData. Drop 'uns' from AnnData to remove this warning, or if you "
                "intend for 'uns' to merge, ensure each AnnData uses unique keys."
            )


@attrs.define(kw_only=True, frozen=True)
class AxisMetadata:
    """Private class"""

    field_name: str  # user-specified join field name
    field_index: pd.Index[Any]  # index of join field values
    enum_values: dict[str, pd.CategoricalDtype]  # dict[col_name, CategoricalDtype]

    @classmethod
    def reduce(cls, ams: list[Self]) -> AxisMetadata:
        assert all(isinstance(a, cls) for a in ams)
        assert all([a.field_name == ams[0].field_name for a in ams])
        if not all(a.field_index.dtype == ams[0].field_index.dtype for a in ams):
            raise SOMAError("All AnnData must have a common dtype for their index.")

        if len(ams) == 0:
            return cls(field_name="", field_index=pd.Index([]), enum_values={})
        if len(ams) == 1:
            return ams[0]

        def _reduce_field_index(indices: list[pd.Index]) -> pd.Index:  # type: ignore[type-arg]
            """reducer for joinid indices"""
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
        enum_values: list[dict[str, pd.CategoricalDtype]],
    ) -> dict[str, pd.CategoricalDtype]:
        """reducer for enum values"""

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
