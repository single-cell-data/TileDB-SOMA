import json
from typing import Any, Dict, Optional, Tuple, Type

import numpy as np
import pyarrow as pa
import shapely
import somacore
from somacore import options

from ._exception import SOMAError


def coordinate_space_from_json(data: str) -> somacore.CoordinateSpace:
    """Returns a coordinate space from a json string."""
    # TODO: Needs good, comprehensive error handling.
    raw = json.loads(data)
    # mypy false positive https://github.com/python/mypy/issues/5313
    return somacore.CoordinateSpace(tuple(somacore.Axis(**axis) for axis in raw))  # type: ignore[misc]


def coordinate_space_to_json(coord_space: somacore.CoordinateSpace) -> str:
    """Returns json string representation of the coordinate space."""
    return json.dumps(
        tuple({"name": axis.name, "unit": axis.unit} for axis in coord_space.axes)
    )


def transform_from_json(data: str) -> somacore.CoordinateTransform:
    """Convert a JSON string representing a CoordinateTransform"""

    raw = json.loads(data)

    try:
        transform_type = raw.pop("transform_type")
    except KeyError:
        raise KeyError(
            "'transform_type' not found when attempting to convert "
            "JSON to CoordinateTransform child class"
        )

    try:
        kwargs = raw.pop("transform")
    except KeyError:
        raise KeyError(
            "'transform' kwargs options not found when attempting to "
            "convert JSON to CoordinateTransform child class"
        )

    coord_transform_init: Dict[str, Type[somacore.CoordinateTransform]] = {
        "AffineTransform": somacore.AffineTransform,
        "ScaleTransform": somacore.ScaleTransform,
        "UniformScaleTransform": somacore.UniformScaleTransform,
        "IdentityTransform": somacore.IdentityTransform,
    }

    try:
        return coord_transform_init[transform_type](**kwargs)
    except KeyError:
        raise KeyError(f"Unrecognized transform type key '{transform_type}'")


def transform_to_json(transform: somacore.CoordinateTransform) -> str:
    """Representing a CoordinateTransform as a JSON string"""

    kwargs: Dict[str, Any] = {
        "input_axes": transform.input_axes,
        "output_axes": transform.output_axes,
    }
    if isinstance(transform, somacore.IdentityTransform):
        pass
    elif isinstance(transform, somacore.UniformScaleTransform):
        kwargs["scale"] = transform.scale
    elif isinstance(transform, somacore.ScaleTransform):
        kwargs["scale_factors"] = transform.scale_factors.tolist()
    elif isinstance(transform, somacore.AffineTransform):
        kwargs["matrix"] = transform.augmented_matrix.tolist()
    else:
        raise TypeError(f"Unrecognized coordinate transform type {type(transform)!r}.")

    transform_type = type(transform).__name__
    return json.dumps({"transform_type": transform_type, "transform": kwargs})


def transform_region(
    region: options.SpatialRegion,
    transform: somacore.CoordinateTransform,
) -> shapely.geometry.base.BaseGeometry:
    if len(transform.input_axes) != 2:
        raise NotImplementedError(
            "Spatial queries are currently only supported for 2D coordinates."
        )
    if isinstance(region, shapely.geometry.base.BaseGeometry):
        if region.has_z:
            raise ValueError("Only 2d shapely geometries are supported.")
        # Following check is currently unneeded, but leaving it for reference if
        # 3D support is added.
        ndim = 3 if region.has_z else 2
        if ndim != len(transform.input_axes):
            raise ValueError(
                f"Input region must have {len(transform.input_axes)} dimension, "
                f"but region with {ndim} dimensions provided."
            )
    else:
        if len(region) != 4:
            raise ValueError(f"Unexpected region with size {len(region)}")
        region = shapely.box(region[0], region[1], region[2], region[3])

    if not isinstance(transform, somacore.AffineTransform):
        raise NotImplementedError("Only affine transforms are supported.")
    aug = transform.augmented_matrix
    affine = aug[:-1, :-1].flatten("C").tolist() + aug[-1, :-1].tolist()
    return shapely.affinity.affine_transform(region, affine)


def process_image_region(
    region: Optional[options.SpatialRegion],
    transform: somacore.CoordinateTransform,
    channel_coords: options.DenseCoord,
    image_type: str,
) -> Tuple[
    options.DenseNDCoords, Optional[options.SpatialRegion], somacore.CoordinateTransform
]:

    if region is None:
        # Select the full region.
        data_region: Optional[options.SpatialRegion] = None
        x_coords: options.DenseCoord = None
        y_coords: options.DenseCoord = None
    else:
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
        x_coords = slice(x_min, x_max)
        y_coords = slice(y_min, y_max)

        # Translate the transform if the region does not start at the origin.
        if x_min != 0 or y_min != 0:
            translate = somacore.AffineTransform(
                transform.output_axes,
                transform.output_axes,
                np.array([[1, 0, -x_min], [0, 1, -y_min], [0, 0, 1]]),
            )

            transform = translate @ transform

    # Get the inverse translation from the data space to the original requested region.
    inv_transform = transform.inverse_transform()

    # Get the dense coordinates for querying the array storing the image.
    coords: options.DenseNDCoords = []
    for axis in image_type:
        if axis == "C":
            coords.append(channel_coords)  # type: ignore[attr-defined]
        if axis == "X":
            coords.append(x_coords)  # type: ignore[attr-defined]
        if axis == "Y":
            coords.append(y_coords)  # type: ignore[attr-defined]
        if axis == "Z":
            raise NotImplementedError(
                "Spatial queries are currently only supported for 2D coordinates."
            )

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
