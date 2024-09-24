import numpy as np
import pytest
from somacore import (
    AffineTransform,
    Axis,
    CoordinateSpace,
    CoordinateTransform,
    IdentityTransform,
)

import tiledbsoma as soma


def check_transform_is_equal(
    actual: CoordinateTransform, desired: CoordinateTransform
) -> None:
    assert actual.input_axes == desired.input_axes
    assert actual.output_axes == desired.output_axes
    if isinstance(desired, IdentityTransform):
        assert isinstance(actual, IdentityTransform)
    elif isinstance(desired, AffineTransform):
        assert isinstance(actual, AffineTransform)
        np.testing.assert_array_equal(actual.augmented_matrix, desired.augmented_matrix)
    else:
        assert False


@pytest.mark.parametrize(
    "original",
    [
        CoordinateSpace((Axis(name="dim0", unit="meter"),)),  # type: ignore[arg-type]
        CoordinateSpace(
            (
                Axis(name="dim0", unit="micrometer"),
                Axis(name="dim1"),
                Axis(name="dim2", unit="micrometer"),
            ),  # type: ignore[arg-type]
        ),
    ],
)
def test_coordinate_system_json_roundtrip(original: CoordinateSpace):
    json_blob = soma._spatial_util.coordinate_space_to_json(original)
    result = soma._spatial_util.coordinate_space_from_json(json_blob)
    assert len(result) == len(original)
    for index in range(len(result)):
        assert result[index] == original[index]


@pytest.mark.parametrize(
    "original",
    [
        IdentityTransform(["y1", "x1"], ["y2", "x2"]),
        AffineTransform(
            ["x1", "y1"],
            ["x2", "y2"],
            np.array([[1.5, 3, 0], [-1.5, 3, 1]], dtype=np.float64),
        ),
    ],
)
def test_transform_json_roundtrip(original: CoordinateTransform):
    json_blob = soma._spatial_util.transform_to_json(original)
    result = soma._spatial_util.transform_from_json(json_blob)
    check_transform_is_equal(result, original)
