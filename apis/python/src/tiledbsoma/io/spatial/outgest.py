# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Any, Dict, Mapping, Sequence, Tuple, Union

import pandas as pd
import somacore

try:
    import spatialdata as sd
except ImportError as err:
    warnings.warn("Experimental spatial outgestor requires the spatialdata package.")
    raise err

try:
    import dask.dataframe as dd
except ImportError as err:
    warnings.warn("Experimental spatial outgestor requires the dask package.")
    raise err

try:
    import geopandas as gpd
except ImportError as err:
    warnings.warn("Experimental spatial outgestor requires the geopandas package.")
    raise err

from ... import Experiment, MultiscaleImage, PointCloudDataFrame, Scene
from ..._constants import SOMA_JOINID, SPATIAL_DISCLAIMER
from ..._spatial_util import transform_from_json
from .. import to_anndata
from ._xarray_backend import dense_nd_array_to_data_array, images_to_datatree

if TYPE_CHECKING:
    from spatialdata.models.models import DataArray, DataTree


def _convert_axis_names(
    coord_axis_names: Tuple[str, ...],
    data_axis_names: Tuple[str, ...] | None = None,
) -> Tuple[Tuple[str, ...], Dict[str, str]]:
    """Convert SOMA axis names to SpatialData axis names.

    Args:
        coord_axis_names: Axis names for the coordinate space.
        data_axis_names: Axis names in the order they are stored for the data.
            If ``None``, set to the same values as ``coord_axis_names``. Defaults
            to ``None``.

    Returns:
        A tuple of axis names to use for SpatialData and a mapping from SOMA axis
        names to SpatialData axis names

    """
    # Create
    ndim = len(coord_axis_names)
    dim_map = {
        "soma_channel": "c",
        coord_axis_names[0]: "x",
        coord_axis_names[1]: "y",
    }
    if ndim >= 3:
        dim_map[coord_axis_names[2]] = "z"
    if data_axis_names is None:
        data_axis_names = coord_axis_names
    spatial_data_axes = tuple(dim_map[axis_name] for axis_name in data_axis_names)
    soma_dim_map = {key: val for key, val in dim_map.items()}
    return spatial_data_axes, soma_dim_map


def _transform_to_spatialdata(
    transform: somacore.CoordinateTransform,
    input_dim_map: Dict[str, str],
    output_dim_map: Dict[str, str],
) -> sd.transformations.BaseTransformation:
    """Returns the equivalent SpatialData transform for a SOMA transform.

    Args:
        transform: The SOMA transform to convert.
        input_dim_map: Mapping from SOMA transform input axes to SpatialData dimension names.
        output_dim_map: Mapping from SOMA transform output axes to SpatialData dimension names.

    Returns:
        Equivalent SpatialData transformation.
    """
    if isinstance(transform, somacore.IdentityTransform):
        return sd.transformations.Identity()
    if isinstance(transform, somacore.ScaleTransform):
        input_axes = tuple(input_dim_map[name] for name in transform.input_axes)
        return sd.transformations.Scale(transform.scale_factors, input_axes)
    if isinstance(transform, somacore.AffineTransform):
        input_axes = tuple(input_dim_map[name] for name in transform.input_axes)
        output_axes = tuple(output_dim_map[name] for name in transform.output_axes)
        return sd.transformations.Affine(
            transform.augmented_matrix, input_axes, output_axes
        )

    raise NotImplementedError(
        f"Support for converting transform of type {type(transform).__name__} is not "
        f"yet implemented."
    )


def to_spatialdata_points(
    points: PointCloudDataFrame,
    *,
    key: str,
    scene_id: str,
    scene_dim_map: Dict[str, str],
    transform: somacore.CoordinateTransform | None,
    soma_joinid_name: str = SOMA_JOINID,
) -> dd.DataFrame:
    """Export a :class:`PointCloudDataFrame` to a :class:`spatialdata.ShapesModel.

    Args:
        points: The point cloud data frame to convert to SpatialData shapes.
        key: Key for the item in the SpatialData object. Used to set a transformation
            to the item itself if no scene transform is provided.
        scene_id: The ID of the scene this point cloud dataframe is from.
        scene_dim_map: A mapping from the axis names of the scene to the corresponding
            SpatialData dimension names.
        transform: The transformation from the coordinate space of the scene this point
            cloud is in to the coordinate space of the point cloud.
        soma_joinid: The name to use for the SOMA joinid.
    """

    # Get the axis names for the spatial data shapes.
    orig_axis_names = points.coordinate_space.axis_names
    new_axis_names, points_dim_map = _convert_axis_names(orig_axis_names)

    # Create the SpatialData transform from the points to the Scene (inverse of the
    # transform SOMA stores).
    if transform is None:
        transforms = {key: sd.transformations.Identity()}
    else:
        transforms = {
            scene_id: _transform_to_spatialdata(
                transform.inverse_transform(), points_dim_map, scene_dim_map
            )
        }

    # Read the pandas dataframe, rename SOMA_JOINID, add metadata, and return.
    df: pd.DataFrame = points.read().concat().to_pandas()
    if soma_joinid_name != SOMA_JOINID:
        df.rename(columns={SOMA_JOINID: soma_joinid_name}, inplace=True)
    return sd.models.PointsModel.parse(df, transformations=transforms)


def to_spatialdata_shapes(
    points: PointCloudDataFrame,
    *,
    key: str,
    scene_id: str,
    scene_dim_map: Dict[str, str],
    transform: somacore.CoordinateTransform | None,
    soma_joinid_name: str = SOMA_JOINID,
) -> gpd.GeoDataFrame:
    """Export a :class:`PointCloudDataFrame` to a :class:`spatialdata.ShapesModel.

    Args:
        points: The point cloud data frame to convert to SpatialData shapes.
        key: Key for the item in the SpatialData object. Used to set a transformation
            to the item itself if no scene transform is provided.
        scene_id: The ID of the scene this point cloud dataframe is from.
        scene_dim_map: A mapping from the axis names of the scene to the corresponding
            SpatialData dimension names.
        transform: The transformation from the coordinate space of the scene this point
            cloud is in to the coordinate space of the point cloud.
        soma_joinid: The name to use for the SOMA joinid.
    """

    # Get the radius for the point cloud.
    try:
        radius = points.metadata["soma_geometry"]
    except KeyError as ke:
        raise KeyError(
            "Missing metadata 'soma_geometry' needed for reading the point cloud "
            "dataframe as a shape."
        ) from ke
    try:
        soma_geometry_type = points.metadata["soma_geometry_type"]
        if soma_geometry_type != "radius":
            raise NotImplementedError(
                f"Support for a point cloud with shape '{soma_geometry_type}' is "
                f"not yet implemented."
            )
    except KeyError as ke:
        raise KeyError("Missing metadata 'soma_geometry_type'.") from ke

    # Get the axis names for the spatial data shapes.
    orig_axis_names = points.coordinate_space.axis_names
    new_axis_names, points_dim_map = _convert_axis_names(orig_axis_names)

    # Create the SpatialData transform from the points to the Scene (inverse of the
    # transform SOMA stores).
    if transform is None:
        transforms = {key: sd.transformations.Identity()}
    else:
        transforms = {
            scene_id: _transform_to_spatialdata(
                transform.inverse_transform(), points_dim_map, scene_dim_map
            )
        }

    data = points.read().concat().to_pandas()
    if soma_joinid_name != SOMA_JOINID:
        data.rename(columns={SOMA_JOINID: soma_joinid_name}, inplace=True)
    data.insert(len(data.columns), "radius", radius)
    ndim = len(orig_axis_names)
    if ndim == 2:
        geometry = gpd.points_from_xy(
            data.pop(orig_axis_names[0]), data.pop(orig_axis_names[1])
        )
    else:
        raise NotImplementedError(
            f"Support for export {ndim}D point cloud dataframes to SpatialData shapes "
            f"is not yet implemented."
        )
    df = gpd.GeoDataFrame(data, geometry=geometry)
    df.attrs["transform"] = transforms
    return df


def to_spatialdata_image(
    image: MultiscaleImage,
    level: Union[str, int] | None = None,
    *,
    key: str,
    scene_id: str,
    scene_dim_map: Dict[str, str],
    transform: somacore.CoordinateTransform | None,
) -> "DataArray":
    """Export a level of a :class:`MultiscaleImage` to a
    :class:`spatialdata.Image2DModel` or :class:`spatialdata.Image3DModel`.

    Args:
        image: The multiscale image to convert to a level from to a SpatialData image.
        level: The level of the multiscale image to convert.
        key: Key for the item in the SpatialData object. Used to set a transformation
            to the item itself if no scene transform is provided.
        scene_id: The ID of the scene this multiscale image is from.
        scene_dim_map: A mapping from the axis names of the scene to the corresponding
            SpatialData dimension names.
        transform: The transformation from the coordinate space of the scene this
            multiscale image is in to the coordinate space of the image itself.
    """
    if not image.has_channel_axis:
        raise NotImplementedError(
            "Support for exporting a MultiscaleImage to without a channel axis to "
            "SpatialData is not yet implemented."
        )

    # Convert from SOMA axis names to SpatialData axis names.
    orig_axis_names = image.coordinate_space.axis_names
    if len(orig_axis_names) not in {2, 3}:
        raise NotImplementedError(
            f"Support for converting a '{len(orig_axis_names)}'D is not yet implemented."
        )
    new_axis_names, image_dim_map = _convert_axis_names(
        orig_axis_names, image.data_axis_order
    )

    # Get the URI of the requested level.
    if level is None:
        if image.level_count != 1:
            raise ValueError(
                "The level must be specified for a multiscale image with more than one "
                "resolution level."
            )
        level = 0
    level_uri = image.level_uri(level)

    if transform is None:
        # Get the transformation from the image level to the highest resolution of the multiscale image.
        scale_transform = image.get_transform_from_level(level)
        sd_transform = _transform_to_spatialdata(
            scale_transform, image_dim_map, image_dim_map
        )
        transformations = {key: sd_transform}
    else:
        # Get the transformation from the image level to the scene:
        # If the result is a single scale transform (or identity transform), output a
        # single transformation. Otherwise, convert to a SpatialData sequence of
        # transformations.
        inv_transform = transform.inverse_transform()
        scale_transform = image.get_transform_from_level(level)
        if isinstance(transform, somacore.ScaleTransform) or isinstance(
            scale_transform, somacore.IdentityTransform
        ):
            # inv_transform @ scale_transform -> applies scale_transform first
            sd_transform = _transform_to_spatialdata(
                inv_transform @ scale_transform, image_dim_map, scene_dim_map
            )
        else:
            sd_transform1 = _transform_to_spatialdata(
                scale_transform, image_dim_map, image_dim_map
            )
            sd_transform2 = _transform_to_spatialdata(
                inv_transform, image_dim_map, scene_dim_map
            )
            # Sequence([sd_transform1, sd_transform2]) -> applies sd_transform1 first
            sd_transform = sd.transformations.Sequence([sd_transform1, sd_transform2])
        transformations = {scene_id: sd_transform}

    # Return array accessor as a dask array.
    return dense_nd_array_to_data_array(
        level_uri,
        dim_names=new_axis_names,
        attrs={"transform": transformations},
        context=image.context,
    )


def to_spatialdata_multiscale_image(
    image: MultiscaleImage,
    *,
    key: str,
    scene_id: str,
    scene_dim_map: Dict[str, str],
    transform: somacore.CoordinateTransform | None,
) -> "DataTree":
    """Export a MultiscaleImage to a DataTree.

    Args:
        image: The multiscale image to convert to a SpatialData image.
        key: Key for the item in the SpatialData object. Used to set a transformation
            to the item itself if no scene transform is provided.
        scene_id: The ID of the scene this multiscale image is from.
        scene_dim_map: A mapping from the axis names of the scene to the corresponding
            SpatialData dimension names.
        transform: The transformation from the coordinate space of the scene this
            multiscale image is in to the coordinate space of the image itself.
    """

    # Check for channel axis.
    if not image.has_channel_axis:
        raise NotImplementedError(
            "Support for exporting a MultiscaleImage to without a channel axis to "
            "SpatialData is not yet implemented."
        )

    # Convert from SOMA axis names to SpatialData axis names.
    orig_axis_names = image.coordinate_space.axis_names
    if len(orig_axis_names) not in {2, 3}:
        raise NotImplementedError(
            f"Support for converting a '{len(orig_axis_names)}'D is not yet implemented."
        )
    new_axis_names, image_dim_map = _convert_axis_names(
        orig_axis_names, image.data_axis_order
    )

    if transform is None:
        spatial_data_transformations = tuple(
            _transform_to_spatialdata(
                image.get_transform_from_level(level),
                image_dim_map,
                image_dim_map,
            )
            for level in range(image.level_count)
        )
    else:
        # Get the transformtion from the image level to the scene:
        # If the result is a single scale transform (or identity transform), output a
        # single transformation. Otherwise, convert to a SpatialData sequence of
        # transformations.
        inv_transform = transform.inverse_transform()
        if isinstance(transform, somacore.ScaleTransform):
            # inv_transform @ scale_transform -> applies scale_transform first
            spatial_data_transformations = tuple(
                _transform_to_spatialdata(
                    inv_transform @ image.get_transform_from_level(level),
                    image_dim_map,
                    scene_dim_map,
                )
                for level in range(image.level_count)
            )
        else:
            sd_scale_transforms = tuple(
                _transform_to_spatialdata(
                    image.get_transform_from_level(level), image_dim_map, image_dim_map
                )
                for level in range(1, image.level_count)
            )
            sd_inv_transform = _transform_to_spatialdata(
                inv_transform, image_dim_map, scene_dim_map
            )

            # First level transform is always the identity, so just directly use
            # inv_transform. For remaining transformations,
            # Sequence([sd_transform1, sd_transform2]) -> applies sd_transform1 first
            spatial_data_transformations = (sd_inv_transform,) + tuple(
                sd.transformations.Sequence([scale_transform, sd_inv_transform])
                for scale_transform in sd_scale_transforms
            )

    # Create a sequence of resolution level.
    image_data_arrays = tuple(
        dense_nd_array_to_data_array(
            uri=image.level_uri(index),
            dim_names=new_axis_names,
            attrs=(
                {"transform": {key: spatial_data_transformations[index]}}
                if transform is None
                else {"transform": {scene_id: spatial_data_transformations[index]}}
            ),
            context=image.context,
        )
        for index, (soma_name, val) in enumerate(image.levels().items())
    )

    return images_to_datatree(image_data_arrays)


def _get_transform_from_collection(
    key: str, metadata: Mapping[str, Any]
) -> somacore.CoordinateTransform | None:
    transform_key = f"soma_scene_registry_{key}"
    if transform_key in metadata:
        transform_json = metadata[transform_key]
        return transform_from_json(transform_json)
    return None


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
    table_kwargs: Mapping[str, Dict[str, Any]] | None = None,
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
