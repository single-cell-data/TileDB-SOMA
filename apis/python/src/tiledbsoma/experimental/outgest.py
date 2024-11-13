# Copyright (c) 2024 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2024 TileDB, Inc
#
# Licensed under the MIT License.
from typing import Dict, Optional, Tuple

import geopandas as gpd
import somacore
import spatialdata as sd

from .. import PointCloudDataFrame
from .._constants import SOMA_JOINID


def _convert_axis_names(
    coord_axis_names: Tuple[str, ...],
    data_axis_names: Optional[Tuple[str, ...]] = None,
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


def _transform_to_spatial_data(
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


def to_spatial_data_shapes(
    points: PointCloudDataFrame,
    *,
    scene_id: str,
    scene_dim_map: Dict[str, str],
    transform: somacore.CoordinateTransform,
    soma_joinid_name: str,
) -> gpd.GeoDataFrame:
    """Export a :class:`PointCloudDataFrame` to a :class:`spatialdata.ShapesModel.

    Args:
        points: The point cloud data frame to convert to SpatialData shapes.
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
    transforms = {
        scene_id: _transform_to_spatial_data(
            transform.inverse_transform(), points_dim_map, scene_dim_map
        )
    }

    data = points.read().concat().to_pandas()
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