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


from ... import Experiment, MultiscaleImage, PointCloudDataFrame, Scene
from ..._constants import SPATIAL_DISCLAIMER
from .. import to_anndata
from ._spatialdata_util import (
    _convert_axis_names,
    _get_transform_from_collection,
    to_spatialdata_image,
    to_spatialdata_multiscale_image,
    to_spatialdata_points,
    to_spatialdata_shapes,
)


def _add_scene_to_spatialdata(
    sdata: sd.SpatialData,
    scene_id: str,
    scene: Scene,
    *,
    obs_id_name: str,
    var_id_name: str,
    measurement_names: Sequence[str] | None = None,
) -> None:
    """Adds items from a Scene to a SpatialData object.

    Args:
        sdata: SpatialData object to update.
        scene_id: Name of the scene that will be added to the SpatialData object.
        scene: The scene that is being added to the SpatialData object.
        obs_id_name: Name to use for the ``soma_joinid`` in ``obsl``.
        var_id_name: Name to use fo the ``soma_joinid in ``varl``.
        measurement_names: The names of measurements to export. If ``None``, all
            measurements are included. Defaults to ``None``.

    """
    # Cannot have spatial data if no coordinate space.
    if scene.coordinate_space is None:
        return

    # Get the map from Scene dimension names to SpatialData dimension names.
    input_axis_names = scene.coordinate_space.axis_names
    _, scene_dim_map = _convert_axis_names(input_axis_names, input_axis_names)

    # Export obsl data to SpatialData.
    if "obsl" in scene:
        for key, df in scene.obsl.items():
            output_key = f"{scene_id}_{key}"
            transform = _get_transform_from_collection(key, scene.obsl.metadata)
            if isinstance(df, PointCloudDataFrame):
                if "soma_geometry" in df.metadata:
                    sdata.shapes[output_key] = to_spatialdata_shapes(
                        df,
                        key=output_key,
                        scene_id=scene_id,
                        scene_dim_map=scene_dim_map,
                        transform=transform,
                        soma_joinid_name=obs_id_name,
                    )
                else:
                    sdata.points[output_key] = to_spatialdata_points(
                        df,
                        key=output_key,
                        scene_id=scene_id,
                        scene_dim_map=scene_dim_map,
                        transform=transform,
                        soma_joinid_name=obs_id_name,
                    )
            else:
                warnings.warn(
                    f"Skipping obsl[{key}] in Scene {scene_id}; unexpected datatype"
                    f" {type(df).__name__}."
                )

    # Export varl data to SpatialData.
    if "varl" in scene:
        for measurement_name, subcoll in scene.varl.items():
            if (
                measurement_names is not None
                and measurement_name not in measurement_names
            ):
                continue
            for key, df in subcoll.items():
                output_key = f"{scene_id}_{measurement_name}_{key}"
                transform = _get_transform_from_collection(key, subcoll.metadata)
                if isinstance(df, PointCloudDataFrame):
                    if "soma_geometry" in df.metadata:
                        sdata.shapes[output_key] = to_spatialdata_shapes(
                            df,
                            key=output_key,
                            scene_id=scene_id,
                            scene_dim_map=scene_dim_map,
                            transform=transform,
                            soma_joinid_name=var_id_name,
                        )
                    else:
                        sdata.points[output_key] = to_spatialdata_points(
                            df,
                            key=output_key,
                            scene_id=scene_id,
                            scene_dim_map=scene_dim_map,
                            transform=transform,
                            soma_joinid_name=var_id_name,
                        )
                else:
                    warnings.warn(
                        f"Skipping varl[{measurement_name}][{key}] in Scene "
                        f"{scene_id}; unexpected datatype {type(df).__name__}."
                    )

    # Export img data to SpatialData.
    if "img" in scene:
        for key, image in scene.img.items():
            output_key = f"{scene_id}_{key}"
            transform = _get_transform_from_collection(key, scene.img.metadata)
            if not isinstance(image, MultiscaleImage):
                warnings.warn(  # type: ignore[unreachable]
                    f"Skipping img[{image}] in Scene {scene_id}; unexpected "
                    f"datatype {type(image).__name__}."
                )
            if image.level_count == 1:
                sdata.images[output_key] = to_spatialdata_image(
                    image,
                    0,
                    key=output_key,
                    scene_id=scene_id,
                    scene_dim_map=scene_dim_map,
                    transform=transform,
                )
            else:
                sdata.images[f"{scene_id}_{key}"] = to_spatialdata_multiscale_image(
                    image,
                    key=output_key,
                    scene_id=scene_id,
                    scene_dim_map=scene_dim_map,
                    transform=transform,
                )


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
    if "ms" in experiment:
        if measurement_names is None:
            ms_keys = tuple(experiment.ms.keys())
        else:
            ms_keys = tuple(measurement_names)
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

    sdata = sd.SpatialData(tables=tables)

    # If no spatial data, return just the tables.
    if "spatial" not in experiment:
        return sdata

    if scene_names is None:
        scene_names = tuple(experiment.spatial.keys())
    else:
        scene_names = tuple(scene_names)

    for scene_id in scene_names:
        scene = experiment.spatial[scene_id]
        _add_scene_to_spatialdata(
            sdata=sdata,
            scene_id=scene_id,
            scene=scene,
            obs_id_name=obs_id_name,
            var_id_name=var_id_name,
            measurement_names=measurement_names,
        )

    return sdata
