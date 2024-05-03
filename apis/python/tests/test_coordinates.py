import numpy as np
import pytest

import tiledbsoma as soma


@pytest.mark.parametrize(
    "original",
    [
        soma.Axis(axis_name="dim0"),
        soma.Axis(axis_name="dim0", axis_type="space"),
        soma.Axis(axis_name="dim0", axis_type="space", axis_unit="meter"),
    ],
)
def test_axis_json_roundtrip(original: soma.Axis):
    json_blob = original.to_json()
    result = soma.Axis.from_json(json_blob)
    assert result.name == original.name
    assert result.type == original.type
    assert result.unit == original.unit


@pytest.mark.parametrize(
    "original",
    [
        soma.CoordinateSystem(
            (soma.Axis(axis_name="dim0", axis_type="space", axis_unit="meter"),)
        ),
        soma.CoordinateSystem(
            (
                soma.Axis(axis_name="dim0", axis_type="space"),
                soma.Axis(axis_name="dim1"),
                soma.Axis(axis_name="dim2", axis_type="space", axis_unit="meter"),
            ),
        ),
        soma.CoordinateSystem([]),
    ],
)
def test_coordinate_system_json_roundtrip(original: soma.CoordinateSystem):
    json_blob = original.to_json()
    result = soma.CoordinateSystem.from_json(json_blob)
    assert len(result) == len(original)
    for index in range(len(result)):
        assert result[index] == original[index]


@pytest.mark.parametrize(
    "original",
    [
        soma.CompositeTransform([soma.ScaleTransform(np.array([1.0, 2.0, 2.0]))]),
    ],
)
def test_coordinate_transform_json_roundtrip(original: soma.CompositeTransform):
    json_blob = original.to_json()
    result = soma.CompositeTransform.from_json(json_blob)
    assert len(result) == len(original)
