import json
from dataclasses import dataclass
from typing import Any, Dict, Optional, Sequence

import anndata as ad
import pandas as pd
from typing_extensions import Self

import tiledbsoma
import tiledbsoma.logging

from .id_mappings import AxisIDMapping, ExperimentIDMapping, get_dataframe_values


@dataclass
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
        return AxisIDMapping(soma_joinids)

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

        if index_field_name is None:
            index_field_name = df.index.name
        if index_field_name is None:
            index_field_name = "index"
        df = df.reset_index()

        data = {}
        next_soma_joinid = 0
        for index in df[index_field_name]:
            data[index] = next_soma_joinid
            next_soma_joinid += 1
        return cls(data, index_field_name)

    def toJSON(self) -> str:
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, indent=4)

    @classmethod
    def fromJSON(cls, s: str) -> Self:
        dikt = json.loads(s)
        return cls(dikt["data"], dikt["field_name"])


@dataclass
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
        """Factory method to compute an label-to-SOMA-join-ID mappings for a single input file in
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
    ) -> Self:
        """Factory method to compute label-to-SOMA-join-ID mappings for a single input file in
        isolation. This is used when a user is ingesting a single AnnData/H5AD to a single SOMA
        experiment, not in append mode, but allowing us to still have the bulk of the ingestor code
        to be non-duplicated between non-append mode and append mode.
        """
        adata = ad.read_h5ad(h5ad_file_name, "r")
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
        # XXX CTX/CFG
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

            with tiledbsoma.Experiment.open(experiment_uri) as exp:
                for batch in exp.obs.read(column_names=["soma_joinid", obs_field_name]):
                    for i in range(len(batch[0])):
                        soma_joinid = batch[0][i].as_py()
                        obs_id = batch[1][i].as_py()
                        if obs_id not in obs_map:
                            obs_map[obs_id] = soma_joinid

                for measurement_name in exp.ms:
                    var_map = {}
                    for batch in exp.ms[measurement_name].var.read(
                        column_names=["soma_joinid", var_field_name]
                    ):
                        for i in range(len(batch[0])):
                            soma_joinid = batch[0][i].as_py()
                            var_id = batch[1][i].as_py()
                            if var_id not in var_map:
                                var_map[var_id] = soma_joinid
                    var_maps[measurement_name] = var_map

            tiledbsoma.logging.logger.info(
                f"Registration: found nobs={len(obs_map)} nvar={len(var_map)} from experiment."
            )

        return cls(
            obs_axis=AxisAmbientLabelMapping(obs_map, obs_field_name),
            var_axes={
                ms_name: AxisAmbientLabelMapping(data, var_field_name)
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
    ) -> Self:
        """Extends registration data to one more AnnData input."""
        tiledbsoma.logging.logger.info("Registration: registering AnnData object.")

        obs_next_soma_joinid = previous.obs_axis.get_next_start_soma_joinid()
        obs_map = previous.obs_axis.data
        obs_ids = get_dataframe_values(adata.obs, obs_field_name)
        for obs_id in obs_ids:
            if obs_id not in obs_map:
                obs_map[obs_id] = obs_next_soma_joinid
                obs_next_soma_joinid += 1

        var_map = previous.var_axes[measurement_name].data
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
                var_maps["raw"] = previous.var_axes["raw"].data

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
            obs_axis=AxisAmbientLabelMapping(obs_map, obs_field_name),
            var_axes={
                ms_name: AxisAmbientLabelMapping(data, var_field_name)
                for ms_name, data in var_maps.items()
            },
        )

    @classmethod
    def from_h5ad_append_on_experiment(
        cls,
        h5ad_file_name: str,
        previous: Self,
        *,
        measurement_name: str,
        obs_field_name: str = "obs_id",
        var_field_name: str = "var_id",
    ) -> Self:
        """Extends registration data to one more H5AD input file."""
        tiledbsoma.logging.logger.info(f"Registration: registering {h5ad_file_name}.")

        adata = ad.read_h5ad(h5ad_file_name, "r")

        return cls.from_anndata_append_on_experiment(
            adata,
            previous,
            measurement_name=measurement_name,
            obs_field_name=obs_field_name,
            var_field_name=var_field_name,
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
    ) -> Self:
        """Extends registration data from the baseline, already-written SOMA
        experiment to include multiple H5AD input files."""

        if experiment_uri is not None and tiledbsoma.Experiment.exists(experiment_uri):
            registration_data = cls.from_isolated_soma_experiment(
                experiment_uri,
                obs_field_name=obs_field_name,
                var_field_name=var_field_name,
            )
        else:
            registration_data = cls(
                obs_axis=AxisAmbientLabelMapping({}, obs_field_name),
                var_axes={
                    measurement_name: AxisAmbientLabelMapping({}, var_field_name),
                    "raw": AxisAmbientLabelMapping({}, var_field_name),
                },
            )

        for h5ad_file_name in h5ad_file_names:
            registration_data = cls.from_h5ad_append_on_experiment(
                h5ad_file_name,
                registration_data,
                measurement_name=measurement_name,
            )

        tiledbsoma.logging.logger.info("Registration: complete.")
        return registration_data

    def show(self) -> None:
        print(f"obs:{len(self.obs_axis.data)}")
        for k, v in self.var_axes.items():
            print(f"{k}/var:{len(v.data)}")

    def toJSON(self) -> str:
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, indent=4)

    @classmethod
    def fromJSON(cls, s: str) -> Self:
        dikt = json.loads(s)
        obs_axis = AxisAmbientLabelMapping(
            dikt["obs_axis"]["data"], dikt["obs_axis"]["field_name"]
        )
        var_axes = {
            k: AxisAmbientLabelMapping(v["data"], v["field_name"])
            for k, v in dikt["var_axes"].items()
        }
        return cls(obs_axis, var_axes)
