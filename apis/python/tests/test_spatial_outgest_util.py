import numpy as np
import pytest
import somacore

soma_outgest = pytest.importorskip("tiledbsoma.io.spatial.outgest")
sd = pytest.importorskip("spatialdata")


@pytest.mark.parametrize(
    "transform, expected",
    [
        (
            somacore.IdentityTransform(("x1", "y1"), ("x2", "y2")),
            sd.transformations.Identity(),
        ),
        (
            somacore.UniformScaleTransform(("x1", "y1"), ("x2", "y2"), 10),
            sd.transformations.Scale([10, 10], ("x", "y")),
        ),
        (
            somacore.ScaleTransform(("x1", "y1"), ("x2", "y2"), [4, 0.1]),
            sd.transformations.Scale([4, 0.1], ("x", "y")),
        ),
        (
            somacore.AffineTransform(
                ["x1", "y1"],
                ["x2", "y2"],
                [[2, 2, 0], [0, 3, 1]],
            ),
            sd.transformations.Affine(
                np.array([[2, 2, 0], [0, 3, 1], [0, 0, 1]]), ("x", "y"), ("x", "y")
            ),
        ),
    ],
)
def test_transform_to_spatialdata(transform, expected):
    input_dim_map = {"x1": "x", "y1": "y", "z1": "z"}
    output_dim_map = {"x2": "x", "y2": "y", "z2": "z"}
    actual = soma_outgest._transform_to_spatialdata(
        transform, input_dim_map, output_dim_map
    )
    assert actual == expected
