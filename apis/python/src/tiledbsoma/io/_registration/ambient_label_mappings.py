import json
from typing import Any, Dict, Optional, Sequence

import anndata as ad
import attrs
import pandas as pd
from typing_extensions import Self

import tiledbsoma
import tiledbsoma.logging
from tiledbsoma.io._util import read_h5ad  # Allow us to read over S3 in backed mode
from tiledbsoma.options import SOMATileDBContext

from .id_mappings import AxisIDMapping, ExperimentIDMapping, get_dataframe_values


@attrs.define(kw_only=True)
class AxisAmbientLabelMapping:
    """
    For all the to-be-appended AnnData/H5AD inputs in SOMA multi-file append-mode ingestion, this
    class tracks the mapping of input-data ``obs`` or ``var`` ID-column name (barcode ID, gene
    symbol) to SOMA join IDs for SOMA experiment ``obs`` or ``var``.

    See module-level comments for more information.
    """

    data: Dict[str, int]
    field_name: str

    def get_next_start_soma_joinid(self) -> int:
        """Once some number of input files have been registered for an ``obs`` or ``var``
        axis, this returned the next as-yet-unused SOMA join ID for the axis."""
        if len(self.data) == 0:
            return 0
        else:
            return max(self.data.values()) + 1

    def id_mapping_from_values(self, input_ids: Sequence[Any]) -> AxisIDMapping:
        """Given registered label-to-SOMA-join-ID mappings for all registered input files for an
        ``obs`` or ``var`` axis, and a list of input-file 0-up offsets, this returns an int-to-int
        mapping from a single input file's ``obs`` or ``var`` axis to the registered SOMA join IDs.
        """
        soma_joinids = []
        for i, input_id in enumerate(input_ids):
            if input_id not in self.data:
                raise ValueError(f"input_id {input_id} not found in registration data")
            soma_joinids.append(self.data[input_id])
        return AxisIDMapping(data=tuple(soma_joinids))

    def id_mapping_from_dataframe(self, df: pd.DataFrame) -> AxisIDMapping:
        """Given registered label-to-SOMA-join-ID mappings for all registered input files for an
        ``obs`` or ``var`` axis, and an input file's dataframe with its 0-up offsets, this returns
        an int-to-int mapping from a single input file's ``obs`` or ``var`` axis to the registered
        SOMA join IDs.
        """
        values = get_dataframe_values(df, self.field_name)
        return self.id_mapping_from_values(values)

    @classmethod
    def from_isolated_dataframe(
        cls,
        df: pd.DataFrame,
        *,
        index_field_name: Optional[str] = None,
    ) -> Self:
        """Factory method to compute an axis label-to-SOMA-join-ID mapping for a single dataframe in
        isolation. This is used when a user is ingesting a single AnnData/H5AD to a single SOMA
        experiment, not in append mode, but allowing us to still have the bulk of the ingestor code
        to be non-duplicated between non-append mode and append mode.
        """
        tiledbsoma.logging.logger.info("Registration: registering AnnData dataframe.")

        index_field_name = index_field_name or df.index.name or "index"
        df = df.reset_index()

        data = {}
        next_soma_joinid = 0
        for index in df[index_field_name]:
            data[index] = next_soma_joinid
            next_soma_joinid += 1
        return cls(data=data, field_name=index_field_name)

    def to_json(self) -> str:
        return json.dumps(self, default=attrs.asdict, sort_keys=True, indent=4)

    @classmethod
    def from_json(cls, s: str) -> Self:
        dikt = json.loads(s)
        return cls(**dikt)


@attrs.define(kw_only=True)
class ExperimentAmbientLabelMapping:
    """
    For all the to-be-appended AnnData/H5AD inputs in SOMA multi-file append-mode ingestion, this
    class contains an ``AxisAmbientLabelMapping`` for ``obs``, and an ``AxisAmbientLabelMapping``
    for ``var`` in each measurement.
    """

    obs_axis: AxisAmbientLabelMapping
    var_axes: Dict[str, AxisAmbientLabelMapping]

    def id_mappings_for_anndata(
        self,
        adata: ad.AnnData,
        *,
        measurement_name: str = "RNA",
        obs_field_name: str = "obs_id",
        var_field_name: str = "var_id",
    ) -> ExperimentIDMapping:
        """
        Given label-to-SOMA-join-ID mappings for all to-be-appended input files, this selects
        out the offset-to-SOMA-join-ID mappings for a single input file.
        """
        obs_axis = self.obs_axis.id_mapping_from_dataframe(adata.obs)
        var_axes = {}
        var_axes[measurement_name] = self.var_axes[
            measurement_name
        ].id_mapping_from_dataframe(adata.var)
        if adata.raw is not None:
            var_axes["raw"] = self.var_axes["raw"].id_mapping_from_dataframe(
                adata.raw.var
            )

        return ExperimentIDMapping(
            obs_axis=obs_axis,
            var_axes=var_axes,
        )

    @classmethod
    def from_isolated_anndata(
        cls,
        adata: ad.AnnData,
        *,
        measurement_name: str,
        obs_field_name: Optional[str] = None,
        var_field_name: Optional[str] = None,
    ) -> Self:
        """Factory method to compute the label-to-SOMA-join-ID mappings for a single input file in
        isolation. This is used when a user is ingesting a single AnnData/H5AD to a single SOMA
        experiment, not in append mode, but allowing us to still have the bulk of the ingestor code
        to be non-duplicated between non-append mode and append mode.
        """
        tiledbsoma.logging.logger.info("Registration: registering AnnData object.")

        obs_axis = AxisAmbientLabelMapping.from_isolated_dataframe(
            adata.obs, index_field_name=obs_field_name
        )
        tiledbsoma.logging.logger.info(
            f"Registration: accumulated to nobs={len(obs_axis.data)}"
        )

        var_axes = {}
        var_axis = AxisAmbientLabelMapping.from_isolated_dataframe(
            adata.var,
            index_field_name=var_field_name,
        )
        var_axes[measurement_name] = var_axis
        tiledbsoma.logging.logger.info(
            f"Registration: accumulated to nvar={len(var_axis.data)}"
        )
        if adata.raw is not None:
            raw_var_axis = AxisAmbientLabelMapping.from_isolated_dataframe(
                adata.raw.var,
                index_field_name=var_field_name,
            )
            var_axes["raw"] = raw_var_axis
            tiledbsoma.logging.logger.info(
                f"Registration: accumulated to nvar={len(raw_var_axis.data)}"
            )

        return cls(
            obs_axis=obs_axis,
            var_axes=var_axes,
        )

    @classmethod
    def from_isolated_h5ad(
        cls,
        h5ad_file_name: str,
        *,
        measurement_name: str,
        obs_field_name: Optional[str] = None,
        var_field_name: Optional[str] = None,
        context: Optional[SOMATileDBContext] = None,
    ) -> Self:
        """Factory method to compute label-to-SOMA-join-ID mappings for a single input file in
        isolation. This is used when a user is ingesting a single AnnData/H5AD to a single SOMA
        experiment, not in append mode, but allowing us to still have the bulk of the ingestor code
        to be non-duplicated between non-append mode and append mode.
        """
        tiledb_ctx = None if context is None else context.tiledb_ctx
        with read_h5ad(h5ad_file_name, mode="r", ctx=tiledb_ctx) as adata:
            return cls.from_isolated_anndata(
                adata,
                measurement_name=measurement_name,
                obs_field_name=obs_field_name,
                var_field_name=var_field_name,
            )

    @classmethod
    def from_isolated_soma_experiment(
        cls,
        experiment_uri: Optional[str] = None,
        *,
        obs_field_name: str = "obs_id",
        var_field_name: str = "var_id",
        context: Optional[SOMATileDBContext] = None,
    ) -> Self:
        """Factory method to compute label-to-SOMA-join-ID mappings for a single SOMA experiment in
        isolation. These are already committed to SOMA storage, so they are the unchangeable inputs
        for a multi-file append-mode ingest to that experiment. The label-to-SOMA-join-ID mappings
        for the input files will be computed on top of this foundation.
        """

        obs_map = {}
        var_maps = {}

        if experiment_uri is None:
            tiledbsoma.logging.logger.info(
                "Registration: starting with empty experiment."
            )

        else:
            tiledbsoma.logging.logger.info(
                f"Registration: starting with experiment {experiment_uri}"
            )

            with tiledbsoma.Experiment.open(experiment_uri, context=context) as exp:
                for batch in exp.obs.read(column_names=["soma_joinid", obs_field_name]):
                    obs_ids = [e.as_py() for e in batch[1]]
                    soma_joinids = [e.as_py() for e in batch[0]]
                    obs_map.update(dict(zip(obs_ids, soma_joinids)))

                for measurement_name in exp.ms:
                    meas = exp.ms[measurement_name]
                    if "var" not in meas:
                        continue
                    expvar = meas.var
                    if var_field_name not in expvar.schema.names:
                        continue
                    var_map = {}
                    for batch in expvar.read(
                        column_names=["soma_joinid", var_field_name]
                    ):
                        var_ids = [e.as_py() for e in batch[1]]
                        soma_joinids = [e.as_py() for e in batch[0]]
                        var_map.update(dict(zip(var_ids, soma_joinids)))
                    var_maps[measurement_name] = var_map

                    tiledbsoma.logging.logger.info(
                        f"Registration: found nobs={len(obs_map)} nvar={len(var_map)} from experiment."
                    )

        return cls(
            obs_axis=AxisAmbientLabelMapping(data=obs_map, field_name=obs_field_name),
            var_axes={
                ms_name: AxisAmbientLabelMapping(data=data, field_name=var_field_name)
                for ms_name, data in var_maps.items()
            },
        )

    @classmethod
    def from_anndata_append_on_experiment(
        cls,
        adata: ad.AnnData,
        previous: Self,
        *,
        measurement_name: str,
        obs_field_name: str = "obs_id",
        var_field_name: str = "var_id",
        append_obsm_varm: bool = False,
    ) -> Self:
        """Extends registration data to one more AnnData input."""
        tiledbsoma.logging.logger.info("Registration: registering AnnData object.")

        # Pre-checks
        if not append_obsm_varm:
            if len(adata.obsm) > 0 or len(adata.varm) > 0:
                raise ValueError(
                    "append-mode ingest of obsm and varm is only supported via explicit opt-in. Please drop them from the inputs, or retry with append_obsm_varm=True."
                )

        if len(adata.obsp) > 0 or len(adata.varp) > 0:
            raise ValueError(
                "append-mode ingest of obsp and varp is not supported. Please retry without them."
            )

        obs_next_soma_joinid = previous.obs_axis.get_next_start_soma_joinid()
        obs_map = dict(previous.obs_axis.data)
        obs_ids = get_dataframe_values(adata.obs, obs_field_name)
        for obs_id in obs_ids:
            if obs_id not in obs_map:
                obs_map[obs_id] = obs_next_soma_joinid
                obs_next_soma_joinid += 1

        var_map = {}
        var_next_soma_joinid = 0
        if measurement_name in previous.var_axes:
            var_map = dict(previous.var_axes[measurement_name].data)
            var_next_soma_joinid = previous.var_axes[
                measurement_name
            ].get_next_start_soma_joinid()
        var_ids = get_dataframe_values(adata.var, var_field_name)
        for var_id in var_ids:
            if var_id not in var_map:
                var_map[var_id] = var_next_soma_joinid
                var_next_soma_joinid += 1

        var_maps = {measurement_name: var_map}

        if adata.raw is None:
            if "raw" in previous.var_axes:
                var_maps["raw"] = dict(previous.var_axes["raw"].data)

        else:
            # One input may not have a raw while the next may have one
            raw_var_map = {}
            raw_var_next_soma_joinid = 0
            if "raw" in previous.var_axes:
                raw_var_axis = previous.var_axes["raw"]
                raw_var_map = raw_var_axis.data
                raw_var_next_soma_joinid = raw_var_axis.get_next_start_soma_joinid()
            raw_var_ids = get_dataframe_values(adata.raw.var, var_field_name)
            for raw_var_id in raw_var_ids:
                if raw_var_id not in raw_var_map:
                    raw_var_map[raw_var_id] = raw_var_next_soma_joinid
                    raw_var_next_soma_joinid += 1

            var_maps["raw"] = raw_var_map

        tiledbsoma.logging.logger.info(
            f"Registration: accumulated to nobs={len(obs_map)} nvar={len(var_map)}."
        )
        return cls(
            obs_axis=AxisAmbientLabelMapping(data=obs_map, field_name=obs_field_name),
            var_axes={
                ms_name: AxisAmbientLabelMapping(data=data, field_name=var_field_name)
                for ms_name, data in var_maps.items()
            },
        )

    @classmethod
    def _acquire_experiment_mappings(
        cls,
        experiment_uri: Optional[str],
        *,
        measurement_name: str,
        obs_field_name: str,
        var_field_name: str,
        context: Optional[SOMATileDBContext] = None,
    ) -> Self:
        """Acquires label-to-ID mappings from the baseline, already-written SOMA experiment."""

        if experiment_uri is not None:
            if not tiledbsoma.Experiment.exists(experiment_uri, context=context):
                raise ValueError("cannot find experiment at URI {experiment_uri}")

            # Pre-check
            with tiledbsoma.Experiment.open(experiment_uri, context=context) as exp:
                if measurement_name not in exp.ms:
                    raise ValueError(
                        f"cannot append: target measurement {measurement_name} is not in experiment {experiment_uri}"
                    )
            registration_data = cls.from_isolated_soma_experiment(
                experiment_uri,
                obs_field_name=obs_field_name,
                var_field_name=var_field_name,
                context=context,
            )
        else:
            registration_data = cls(
                obs_axis=AxisAmbientLabelMapping(data={}, field_name=obs_field_name),
                var_axes={
                    measurement_name: AxisAmbientLabelMapping(
                        data={}, field_name=var_field_name
                    ),
                    "raw": AxisAmbientLabelMapping(data={}, field_name=var_field_name),
                },
            )
        return registration_data

    @classmethod
    def from_anndata_appends_on_experiment(
        cls,
        experiment_uri: Optional[str],
        adatas: Sequence[ad.AnnData],
        *,
        measurement_name: str,
        obs_field_name: str,
        var_field_name: str,
        append_obsm_varm: bool = False,
        context: Optional[SOMATileDBContext] = None,
    ) -> Self:
        """Extends registration data from the baseline, already-written SOMA
        experiment to include multiple H5AD input files. If ``experiment_uri``
        is ``None`` then you will be computing registrations only for the input
        ``AnnData`` objects. If ``experiment_uri`` is not ``None`` then it is
        an error if the experiment is not accessible."""

        registration_data = cls._acquire_experiment_mappings(
            experiment_uri,
            measurement_name=measurement_name,
            obs_field_name=obs_field_name,
            var_field_name=var_field_name,
            context=context,
        )

        for adata in adatas:
            registration_data = cls.from_anndata_append_on_experiment(
                adata,
                registration_data,
                measurement_name=measurement_name,
                append_obsm_varm=append_obsm_varm,
                obs_field_name=obs_field_name,
                var_field_name=var_field_name,
            )

        tiledbsoma.logging.logger.info("Registration: complete.")
        return registration_data

    @classmethod
    def from_h5ad_append_on_experiment(
        cls,
        h5ad_file_name: str,
        previous: Self,
        *,
        measurement_name: str,
        obs_field_name: str = "obs_id",
        var_field_name: str = "var_id",
        append_obsm_varm: bool = False,
        context: Optional[SOMATileDBContext] = None,
    ) -> Self:
        """Extends registration data to one more H5AD input file."""
        tiledbsoma.logging.logger.info(f"Registration: registering {h5ad_file_name}.")

        tiledb_ctx = None if context is None else context.tiledb_ctx
        with read_h5ad(h5ad_file_name, mode="r", ctx=tiledb_ctx) as adata:
            return cls.from_anndata_append_on_experiment(
                adata,
                previous,
                measurement_name=measurement_name,
                obs_field_name=obs_field_name,
                var_field_name=var_field_name,
                append_obsm_varm=append_obsm_varm,
            )

    @classmethod
    def from_h5ad_appends_on_experiment(
        cls,
        experiment_uri: Optional[str],
        h5ad_file_names: Sequence[str],
        *,
        measurement_name: str,
        obs_field_name: str,
        var_field_name: str,
        append_obsm_varm: bool = False,
        context: Optional[SOMATileDBContext] = None,
    ) -> Self:
        """Extends registration data from the baseline, already-written SOMA
        experiment to include multiple H5AD input files."""

        registration_data = cls._acquire_experiment_mappings(
            experiment_uri,
            measurement_name=measurement_name,
            obs_field_name=obs_field_name,
            var_field_name=var_field_name,
            context=context,
        )

        for h5ad_file_name in h5ad_file_names:
            registration_data = cls.from_h5ad_append_on_experiment(
                h5ad_file_name,
                registration_data,
                measurement_name=measurement_name,
                append_obsm_varm=append_obsm_varm,
                obs_field_name=obs_field_name,
                var_field_name=var_field_name,
                context=context,
            )

        tiledbsoma.logging.logger.info("Registration: complete.")
        return registration_data

    def __str__(self) -> str:
        lines = [f"obs:{len(self.obs_axis.data)}"]
        for k, v in self.var_axes.items():
            lines.append(f"{k}/var:{len(v.data)}")
        return "\n".join(lines)

    def to_json(self) -> str:
        return json.dumps(self, default=attrs.asdict, sort_keys=True, indent=4)

    @classmethod
    def from_json(cls, s: str) -> Self:
        dikt = json.loads(s)
        obs_axis = AxisAmbientLabelMapping(
            data=dikt["obs_axis"]["data"], field_name=dikt["obs_axis"]["field_name"]
        )
        var_axes = {
            k: AxisAmbientLabelMapping(data=v["data"], field_name=v["field_name"])
            for k, v in dikt["var_axes"].items()
        }
        return cls(obs_axis=obs_axis, var_axes=var_axes)
