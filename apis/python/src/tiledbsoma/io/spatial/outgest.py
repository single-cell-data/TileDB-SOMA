# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

from __future__ import annotations

import warnings
from typing import Any, Mapping, Sequence

try:
    import spatialdata as sd
except ImportError as err:
    warnings.warn("Experimental spatial outgestor requires the spatialdata package.")
    raise err


from ... import Experiment
from ..._constants import SPATIAL_DISCLAIMER
from .. import to_anndata
from ._spatialdata_util import _spatial_to_spatialdata


def to_spatialdata(
    experiment: Experiment,
    *,
    measurement_names: Sequence[str] | None = None,
    scene_names: Sequence[str] | None = None,
    obs_id_name: str = "obs_id",
    var_id_name: str = "var_id",
    table_kwargs: Mapping[str, dict[str, Any]] | None = None,
) -> sd.SpatialData:
    """Converts the experiment group to SpatialData format.

    Args:
        experiment:
        measurement_names: The names of measurements to export. If ``None``, all
            measurements are included. Defaults to ``None``.
        scene_names: The names of the scenes to export. If ``None``, all scenes are
            included. Defaults to ``None``.
        obs_id_name: Column name to use for ``obs`` dataframes. Defaults to
            ``"obs_id"``.
        var_id_name: Column name to use for ``var`` dataframes. Defaults to
            ``"var_id|``.
        table_kwargs: Optional mapping from measurment name to keyword arguments to
            pass to table conversions. See :method:`to_anndata` for possible keyword
            arguments.
    """
    warnings.warn(SPATIAL_DISCLAIMER)

    # Read non-spatial data into Anndata tables.
    if measurement_names is not None:
        measurement_names = tuple(measurement_names)
    if "ms" in experiment:
        if measurement_names is None:
            ms_keys = tuple(experiment.ms.keys())
        else:
            ms_keys = measurement_names
        if table_kwargs is None:
            table_kwargs = {}
        tables = {
            measurement_name: (
                to_anndata(
                    experiment,
                    measurement_name,
                    obs_id_name=obs_id_name,
                    var_id_name=var_id_name,
                    **table_kwargs[measurement_name],
                )
                if measurement_name in table_kwargs
                else to_anndata(
                    experiment,
                    measurement_name,
                    obs_id_name=obs_id_name,
                    var_id_name=var_id_name,
                )
            )
            for measurement_name in ms_keys
        }
    else:
        tables = {}

    # If no spatial data, return just the tables.
    if "spatial" not in experiment:
        return sd.SpatialData(tables=tables)

    if scene_names is None:
        scene_names = tuple(experiment.spatial.keys())
    else:
        scene_names = tuple(scene_names)

    return _spatial_to_spatialdata(
        experiment.spatial,
        tables,
        scene_names=scene_names,
        measurement_names=measurement_names,
        obs_id_name=obs_id_name,
        var_id_name=var_id_name,
    )
