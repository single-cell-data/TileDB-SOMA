import json

from somacore import (
    AffineTransform,
    Axis,
    CoordinateSpace,
    CoordinateTransform,
    IdentityTransform,
)


def coordinate_space_from_json(data: str) -> CoordinateSpace:
    """Returns a coordinate space from a json string."""
    # TODO: Needs good, comprehensive error handling.
    raw = json.loads(data)
    return CoordinateSpace(tuple(Axis(**axis) for axis in raw))


def coordinate_space_to_json(coord_space: CoordinateSpace) -> str:
    """Returns json string representation of the coordinate space."""
    return json.dumps(
        tuple(
            {"name": axis.name, "units": axis.units, "scale": axis.scale}
            for axis in coord_space.axes
        )
    )


def transform_from_json(data: str) -> CoordinateTransform:
    """TODO: Add docstring"""
    raw = json.loads(data)
    try:
        transform_type = raw.pop("transform_type")
    except KeyError:
        raise KeyError()  # TODO Add error message
    try:
        kwargs = raw.pop("transform")
    except KeyError:
        raise KeyError()  # TODO Add error message
    if transform_type == "IdentityTransform":
        return IdentityTransform(**kwargs)
    elif transform_type == "AffineTransform":
        return AffineTransform(**kwargs)
    else:
        raise KeyError("Unrecognized transform type key 'transform_type'")


def transform_to_json(transform: CoordinateTransform) -> str:
    kwargs = {
        "input_axes": transform.input_axes,
        "output_axes": transform.output_axes,
    }
    if isinstance(transform, IdentityTransform):
        pass
    elif isinstance(transform, AffineTransform):
        kwargs["matrix"] = transform.augmented_matrix.tolist()
    else:
        raise TypeError(f"Unrecognized coordinate transform type {type(transform)!r}.")

    transform_type = type(transform).__name__
    return json.dumps({"transform_type": transform_type, "transform": kwargs})
