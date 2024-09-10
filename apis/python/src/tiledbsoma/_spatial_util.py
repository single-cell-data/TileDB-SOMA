from typing import Tuple

import numpy as np
import shapely
import somacore
from somacore import options


def transform_region(
    region: options.SpatialRegion,
    transform: somacore.CoordinateTransform,
) -> options.SpatialRegion:
    if transform.input_rank != 2:
        raise NotImplementedError(
            "Spatial queries are currently only supported for 2D coordinates."
        )
    if isinstance(region, shapely.GeometryType):
        if region.has_z:
            raise ValueError("Only 2d shapely geometries are supported.")
        ndim = 3 if region.has_z else 2
        if ndim != transform.input_rank:
            raise ValueError(
                "Input region must have {len(transform.input_rank)} dimension, but "
                "region with {ndim} dimensions provided."
            )
        if not isinstance(transform, somacore.AffineTransform):
            raise NotImplementedError("Only affine transforms are supported.")
        aug = transform.augmented_matrix
        affine = aug[:-1, :-1].flatten("C").tolist() + aug[-1, :-1].tolist()
        return shapely.affine_transform(affine)
    region = np.array(region)
    ndim, ncoords = region.shape
    if region.ndim != 2:
        raise ValueError(f"Provided region has unexpected shape={region.shape}.")
    if ncoords != 4:
        raise NotImplementedError("Numpy array describe a box with shape {(2, 4)}.")
    if len(region) != transform.input_rank:
        raise ValueError(
            "Input region must have {len(transform.input_rank)} dimension, but "
            "{len(region) dimensions provided."
        )
    # TODO: Update to work for any affine transforms.
    if not isinstance(transform, somacore.ScaleTransform):
        raise NotImplementedError("Numpy region only supported for scale transforms.")
    return transform.augmented_matrix[:-1, :-1] @ region


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
    if isinstance(region, shapely.GeometryType):
        (x_min, y_min, x_max, y_max) = region
    else:
        (x_min, y_min, x_max, y_max) = np.min(
            region[0, :], np.min(region[1, :]), np.max(region[0, :]), np.max(region[1,])
        )
    x_min = max(0, int(np.floor(x_min)))
    y_min = max(0, int(np.floor(y_min)))
    x_max = int(np.ceil(x_max))
    y_max = int(np.ceil(y_max))

    # Get the inverse translation from the data space to the original requested region.
    if x_min != 0 or y_min != 0:
        translate = somacore.AffineTransform(
            transform.output_axes,
            transform.output_axes,
            np.array([[1, 0, x_min], [0, 1, y_min], [0, 0, 1]]),
        )
        transform = transform * translate
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
