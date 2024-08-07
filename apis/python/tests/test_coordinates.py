import numpy as np
import pytest

import tiledbsoma as soma


@pytest.mark.parametrize(
    "original",
    [
        soma.Axis(name="dim0"),
        soma.Axis(name="dim0", units="micrometer"),
        soma.Axis(name="dim0", units="nanometer", scale=np.float64(65.0)),
    ],
)
def test_axis_json_roundtrip(original: soma.Axis):
    json_blob = original.to_json()
    result = soma.Axis.from_json(json_blob)
    assert result.name == original.name
    assert result.units == original.units
    assert result.scale == original.scale


@pytest.mark.parametrize(
    "original",
    [
        soma.CoordinateSpace((soma.Axis(name="dim0", units="meter"),)),
        soma.CoordinateSpace(
            (
                soma.Axis(name="dim0", units="micrometer"),
                soma.Axis(name="dim1"),
                soma.Axis(name="dim2", units="micrometer", scale=np.float64(65.0)),
            ),
        ),
    ],
)
def test_coordinate_system_json_roundtrip(original: soma.CoordinateSpace):
    json_blob = original.to_json()
    result = soma.CoordinateSpace.from_json(json_blob)
    assert len(result) == len(original)
    for index in range(len(result)):
        assert result[index] == original[index]


@pytest.mark.parametrize(
    "original",
    [
        soma.IdentityCoordinateTransform(["y1", "x1"], ["y2", "x2"]),
        soma.AffineCoordinateTransform(
            ["x", "y"],
            ["a", "b", "c"],
            np.array([[1.5, 3.0], [-1.5, 3.0], [1.0, 0]], dtype=np.float64),
        ),
    ],
)
def test_transform_json_roundtrip(original: soma._coordinates.CoordinateTransform):
    json_blob = soma._coordinates.transform_to_json(original)
    result = soma._coordinates.transform_from_json(json_blob)
    assert original.input_axes == result.input_axes
    assert original.output_axes == result.output_axes
    if isinstance(original, soma.IdentityCoordinateTransform):
        assert isinstance(result, soma.IdentityCoordinateTransform)
    elif isinstance(original, soma.AffineCoordinateTransform):
        assert isinstance(result, soma.AffineCoordinateTransform)
        np.testing.assert_array_equal(result._affine_matrix, original._affine_matrix)
    else:
        assert False
