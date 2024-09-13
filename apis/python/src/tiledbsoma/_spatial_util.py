from typing import Dict, Optional, Tuple

import numpy as np
import pyarrow as pa
import shapely
import somacore
from somacore import options

from ._exception import SOMAError


def transform_region(
    region: options.SpatialRegion,
    transform: somacore.CoordinateTransform,
) -> shapely.GeometryType:
    if transform.input_rank != 2:
        raise NotImplementedError(
            "Spatial queries are currently only supported for 2D coordinates."
        )
    if isinstance(region, shapely.GeometryType):
        if region.has_z:
            raise ValueError("Only 2d shapely geometries are supported.")
        # Following check is currently unneeded, but leaving it for reference if
        # 3D support is added.
        ndim = 3 if region.has_z else 2
        if ndim != transform.input_rank:
            raise ValueError(
                "Input region must have {len(transform.input_rank)} dimension, but "
                "region with {ndim} dimensions provided."
            )
    else:
        if len(region) != 4:
            raise ValueError("Unexpected region with size {len(region)}")
        region = shapely.box(region[0], region[1], region[2], region[3])

    if not isinstance(transform, somacore.AffineTransform):
        raise NotImplementedError("Only affine transforms are supported.")
    aug = transform.augmented_matrix
    affine = aug[:-1, :-1].flatten("C").tolist() + aug[-1, :-1].tolist()
    return shapely.affinity.affine_transform(region, affine)


def process_image_region(
    region: options.SpatialRegion,
    transform: somacore.CoordinateTransform,
    channel_coords: options.DenseCoord,
    image_type: str,
) -> Tuple[options.DenseNDCoords, options.SpatialRegion, somacore.CoordinateTransform]:
    # Get the transformed region the user is selecting in the data space.
    # Note: transform region verifies only 2D data. This function is hard-coded to
    # make the same assumption.
    data_region = transform_region(region, transform)

    # Convert the region to a bounding box. Round values of bounding box to integer
    # values. Include any partially intersected pixels.
    (x_min, y_min, x_max, y_max) = shapely.bounds(data_region)
    x_min = max(0, int(np.floor(x_min)))
    y_min = max(0, int(np.floor(y_min)))
    x_max = int(np.ceil(x_max))
    y_max = int(np.ceil(y_max))

    # Get the inverse translation from the data space to the original requested region.
    if x_min != 0 or y_min != 0:
        translate = somacore.AffineTransform(
            transform.output_axes,
            transform.output_axes,
            np.array([[1, 0, -x_min], [0, 1, -y_min], [0, 0, 1]]),
        )
        transform = translate @ transform
    inv_transform = transform.inverse_transform()

    coords: options.DenseNDCoords = []
    for axis in image_type:
        if axis == "C":
            coords.append(channel_coords)  # type: ignore[attr-defined]
        if axis == "X":
            coords.append(slice(x_min, x_max))  # type: ignore[attr-defined]
        if axis == "Y":
            coords.append(slice(y_min, y_max))  # type: ignore[attr-defined]
        if axis == "Z":
            raise NotImplementedError()
    return (coords, data_region, inv_transform)


def process_spatial_df_region(
    region: Optional[options.SpatialRegion],
    transform: somacore.CoordinateTransform,
    coords_by_name: Dict[str, options.SparseDFCoord],
    index_columns: Tuple[str, ...],
    axis_names: Tuple[str, ...],
    schema: pa.Schema,
) -> Tuple[
    options.SparseDFCoords,
    Optional[options.SpatialRegion],
    somacore.CoordinateTransform,
]:
    # Check provided coords are valid.
    if not set(axis_names).isdisjoint(coords_by_name):
        raise KeyError("Extra coords cannot contain a spatial index column.")
    if not set(index_columns).issuperset(coords_by_name):
        raise KeyError("Extra coords must be index columns.")

    # Transform the region into the data region and add the spatial coordinates
    # to the coords_by_name map.
    if region is None:
        # Leave spatial coords as None - this will select the entire region.
        data_region: Optional[options.SpatialRegion] = None
    else:
        # Restricted to guarantee data region is a box.
        if isinstance(region, shapely.GeometryType):
            raise NotImplementedError(
                "Support for querying point clouds by geometries is not yet implemented."
            )
        if not isinstance(transform, somacore.ScaleTransform):
            raise NotImplementedError(
                f"Support for querying point clouds with a transform of type "
                f"{type(transform)!r} our a bounding box region is not yet supported."
            )
        # Note: transform_region currently only supports 2D regions. This code block
        # operates under the same assumption.
        data_region = transform_region(region, transform)

        def axis_slice(a: float, b: float, dtype: pa.DataType) -> slice:
            if pa.types.is_floating(dtype):
                return slice(a, b)
            if pa.types.is_integer(dtype):
                # Round so only points within the requested region are returned.
                return slice(int(np.ceil(a)), int(np.floor(b)))
            raise SOMAError(f"Unexpected spatial axis datatype {dtype}")

        # Add the transform region to coords. This gets the bounding box for the
        # requested region.
        (x_min, y_min, x_max, y_max) = shapely.bounds(data_region)
        coords_by_name[axis_names[0]] = axis_slice(
            x_min, x_max, schema.field(axis_names[0]).type
        )
        coords_by_name[axis_names[1]] = axis_slice(
            y_min, y_max, schema.field(axis_names[1]).type
        )

    coords = tuple(coords_by_name.get(index_name) for index_name in index_columns)
    inv_transform = transform.inverse_transform()
    return (coords, data_region, inv_transform)
