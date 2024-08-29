import numpy as np
import pytest
from somacore import Axis, CoordinateSpace

import tiledbsoma as soma


def check_transform_is_equal(
    actual: soma.CoordinateTransform, desired: soma.CoordinateTransform
) -> None:
    assert actual.input_axes == desired.input_axes
    assert actual.output_axes == desired.output_axes
    if isinstance(desired, soma.IdentityCoordinateTransform):
        assert isinstance(actual, soma.IdentityCoordinateTransform)
    elif isinstance(desired, soma.AffineCoordinateTransform):
        assert isinstance(actual, soma.AffineCoordinateTransform)
        np.testing.assert_array_equal(actual.augmented_matrix, desired.augmented_matrix)
    else:
        assert False


@pytest.mark.parametrize(
    "original",
    [
        CoordinateSpace((Axis(name="dim0", units="meter"),)),
        CoordinateSpace(
            (
                Axis(name="dim0", units="micrometer"),
                Axis(name="dim1"),
                Axis(name="dim2", units="micrometer", scale=np.float64(65.0)),
            ),
        ),
    ],
)
def test_coordinate_system_json_roundtrip(original: CoordinateSpace):
    json_blob = soma._coordinates.coordinate_space_to_json(original)
    result = soma._coordinates.coordinate_space_from_json(json_blob)
    assert len(result) == len(original)
    for index in range(len(result)):
        assert result[index] == original[index]


@pytest.mark.parametrize(
    "original",
    [
        soma.IdentityCoordinateTransform(["y1", "x1"], ["y2", "x2"]),
        soma.AffineCoordinateTransform(
            ["x1", "y1"],
            ["x2", "y2"],
            np.array([[1.5, 3, 0], [-1.5, 3, 1]], dtype=np.float64),
        ),
    ],
)
def test_transform_json_roundtrip(original: soma._coordinates.CoordinateTransform):
    json_blob = soma._coordinates.transform_to_json(original)
    result = soma._coordinates.transform_from_json(json_blob)
    check_transform_is_equal(result, original)


@pytest.mark.parametrize(
    ("input", "expected"),
    [
        (
            soma.AffineCoordinateTransform(
                ["x1", "y1"],
                ["x2", "y2"],
                [[2, 2, 0], [0, 3, 1]],
            ),
            np.array([[2, 2, 0], [0, 3, 1], [0, 0, 1]], np.float64),
        )
    ],
)
def test_affine_augmented_matrix(input, expected):
    result = input.augmented_matrix
    np.testing.assert_array_equal(result, expected)


@pytest.mark.parametrize(
    ("transform_a", "transform_b", "expected"),
    [
        (
            soma.IdentityCoordinateTransform(["x1", "y1"], ["x2", "y2"]),
            soma.IdentityCoordinateTransform(["x2", "y2"], ["x3", "y3"]),
            soma.IdentityCoordinateTransform(["x1", "y1"], ["x3", "y3"]),
        ),
        (
            soma.IdentityCoordinateTransform(["x1", "y1"], ["x2", "y2"]),
            soma.AffineCoordinateTransform(
                ["x2", "y2"],
                ["x3", "y3"],
                np.array([[1.5, 3.0, 0.0], [-1.5, 3.0, 1.0]], dtype=np.float64),
            ),
            soma.AffineCoordinateTransform(
                ["x1", "y1"],
                ["x3", "y3"],
                np.array([[1.5, 3.0, 0.0], [-1.5, 3.0, 1.0]], dtype=np.float64),
            ),
        ),
        (
            soma.AffineCoordinateTransform(
                ["x1", "y1"],
                ["x2", "y2"],
                np.array([[1.5, 3.0, 0.0], [-1.5, 3.0, 1.0]], dtype=np.float64),
            ),
            soma.IdentityCoordinateTransform(["x2", "y2"], ["x3", "y3"]),
            soma.AffineCoordinateTransform(
                ["x1", "y1"],
                ["x3", "y3"],
                np.array([[1.5, 3.0, 0.0], [-1.5, 3.0, 1.0]], dtype=np.float64),
            ),
        ),
        (
            soma.AffineCoordinateTransform(
                ["x1", "y1"],
                ["x2", "y2"],
                np.array([[2.0, 0.0, 1.0], [0.0, 4.0, 1.0]], dtype=np.float64),
            ),
            soma.AffineCoordinateTransform(
                ["x2", "y2"],
                ["x3", "y3"],
                np.array([[1.0, 1.0, -1.0], [0.0, 1.0, 2.0]], dtype=np.float64),
            ),
            soma.AffineCoordinateTransform(
                ["x1", "y1"],
                ["x3", "y3"],
                np.array([[2.0, 2.0, -1.0], [0.0, 4.0, 9.0]], dtype=np.float64),
            ),
        ),
    ],
)
def test_multiply_tranform(
    transform_a,
    transform_b,
    expected: soma.CoordinateTransform,
):
    result = transform_a * transform_b
    check_transform_is_equal(result, expected)
