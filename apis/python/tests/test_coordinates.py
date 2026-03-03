import numpy as np
import pytest

from tiledbsoma import (
    AffineTransform,
    Axis,
    CoordinateSpace,
    CoordinateTransform,
    IdentityTransform,
    ScaleTransform,
    UniformScaleTransform,
)


def check_transform_is_equal(actual: CoordinateTransform, desired: CoordinateTransform) -> None:
    assert actual.input_axes == desired.input_axes
    assert actual.output_axes == desired.output_axes
    if isinstance(desired, IdentityTransform):
        assert isinstance(actual, IdentityTransform)
    elif isinstance(desired, UniformScaleTransform):
        assert isinstance(actual, UniformScaleTransform)
        assert actual.scale == desired.scale
    elif isinstance(desired, ScaleTransform):
        assert isinstance(actual, ScaleTransform)
        np.testing.assert_array_equal(actual.scale_factors, desired.scale_factors)
    elif isinstance(desired, AffineTransform):
        assert isinstance(actual, AffineTransform)
        np.testing.assert_array_equal(actual.augmented_matrix, desired.augmented_matrix)
    else:
        assert False


def test_invalid_axis_name():
    with pytest.raises(ValueError):
        Axis("soma_axis")


def test_coordinate_space():
    coord_space = CoordinateSpace(
        (Axis("x", unit="nanometer"), Axis("y", unit="nanometer"))  # type: ignore[arg-type]
    )
    assert len(coord_space) == 2
    assert coord_space.axis_names == ("x", "y")
    assert coord_space[0] == Axis("x", unit="nanometer")


def test_coordiante_space_from_axis_names():
    coord_space = CoordinateSpace.from_axis_names(["alpha", "beta"])
    assert len(coord_space) == 2
    assert coord_space.axis_names == ("alpha", "beta")
    assert coord_space[0] == Axis("alpha", unit=None)
    assert coord_space[1] == Axis("beta", unit=None)


@pytest.mark.parametrize(
    ("input", "expected"),
    [
        (
            AffineTransform(
                ["x1", "y1"],
                ["x2", "y2"],
                [[2, 2, 0], [0, 3, 1]],
            ),
            np.array([[2, 2, 0], [0, 3, 1], [0, 0, 1]], np.float64),
        ),
        (
            AffineTransform(
                ["x1", "y1"],
                ["x2", "y2"],
                [[2, 2], [0, 3]],
            ),
            np.array([[2, 2, 0], [0, 3, 0], [0, 0, 1]], np.float64),
        ),
        (
            AffineTransform(
                ["x1", "y1"],
                ["x2", "y2"],
                [[2, 2, 0], [0, 3, 1], [0, 0, 1]],
            ),
            np.array([[2, 2, 0], [0, 3, 1], [0, 0, 1]], np.float64),
        ),
    ],
)
def test_affine_augmented_matrix(input, expected):
    result = input.augmented_matrix
    np.testing.assert_array_equal(result, expected)


@pytest.mark.parametrize(("input_matrix",), [([1, 2, 3],), ([[1, 0, 1], [0, 1, 1], [1, 0, 1]],)])
def test_affine_matrix_value_error(input_matrix):
    with pytest.raises(ValueError):
        AffineTransform(("x1", "y1"), ("x2", "y2"), input_matrix)


def test_bad_number_of_scale_factors():
    with pytest.raises(ValueError):
        ScaleTransform(("x1", "y1"), ("x2", "y2"), [1, 2, 3])


@pytest.mark.parametrize(
    ("input", "expected"),
    [
        (
            AffineTransform(
                ["x1", "y1"],
                ["x2", "y2"],
                [[1, 0, 0], [0, 1, 0]],
            ),
            AffineTransform(
                ["x2", "y2"],
                ["x1", "y1"],
                [[1, 0, 0], [0, 1, 0]],
            ),
        ),
        (
            AffineTransform(
                ["x1", "y1"],
                ["x2", "y2"],
                [[1, 0, 5], [0, 1, 10]],
            ),
            AffineTransform(
                ["x2", "y2"],
                ["x1", "y1"],
                [[1, 0, -5], [0, 1, -10]],
            ),
        ),
        (
            AffineTransform(
                ["x1", "y1"],
                ["x2", "y2"],
                [[2, 0, -5], [0, 4, 5]],
            ),
            AffineTransform(
                ["x2", "y2"],
                ["x1", "y1"],
                [[0.5, 0, 2.5], [0, 0.25, -1.25]],
            ),
        ),
        (
            ScaleTransform(["x1", "y1"], ["x2", "y2"], [4, 0.1]),
            ScaleTransform(["x2", "y2"], ["x1", "y1"], [0.25, 10]),
        ),
        (
            UniformScaleTransform(["x1", "y1"], ["x2", "y2"], 10),
            UniformScaleTransform(["x2", "y2"], ["x1", "y1"], 0.1),
        ),
        (
            IdentityTransform(["x1", "y1"], ["x2", "y2"]),
            IdentityTransform(["x2", "y2"], ["x1", "y1"]),
        ),
    ],
)
def test_inverse_transform(input, expected):
    result = input.inverse_transform()
    check_transform_is_equal(result, expected)
    result_matrix = input.augmented_matrix @ result.augmented_matrix
    expected_matrix: np.ndarray = np.identity(len(result.input_axes) + 1, dtype=np.float64)
    np.testing.assert_allclose(result_matrix, expected_matrix)


def test_uniform_scale_factor():
    UniformScaleTransform(["x1", "y1"], ["x2", "y2"], 1.5)
    UniformScaleTransform(["x1", "y1"], ["x3", "y3"], 1.5)


@pytest.mark.parametrize(
    ("transform_a", "transform_b", "expected"),
    [
        (
            IdentityTransform(["x2", "y2"], ["x3", "y3"]),
            IdentityTransform(["x1", "y1"], ["x2", "y2"]),
            IdentityTransform(["x1", "y1"], ["x3", "y3"]),
        ),
        (
            IdentityTransform(["x2", "y2"], ["x3", "y3"]),
            UniformScaleTransform(["x1", "y1"], ["x2", "y2"], 1.5),
            UniformScaleTransform(["x1", "y1"], ["x3", "y3"], 1.5),
        ),
        (
            IdentityTransform(["x2", "y2"], ["x3", "y3"]),
            ScaleTransform(["x1", "y1"], ["x2", "y2"], np.array([1.5, 3.0], dtype=np.float64)),
            ScaleTransform(["x1", "y1"], ["x3", "y3"], np.array([1.5, 3.0], dtype=np.float64)),
        ),
        (
            IdentityTransform(["x2", "y2"], ["x3", "y3"]),
            AffineTransform(
                ["x1", "y1"],
                ["x2", "y2"],
                np.array([[1.5, 3.0, 0.0], [-1.5, 3.0, 1.0]], dtype=np.float64),
            ),
            AffineTransform(
                ["x1", "y1"],
                ["x3", "y3"],
                np.array([[1.5, 3.0, 0.0], [-1.5, 3.0, 1.0]], dtype=np.float64),
            ),
        ),
        (
            UniformScaleTransform(["x2", "y2"], ["x3", "y3"], 1.5),
            IdentityTransform(["x1", "y1"], ["x2", "y2"]),
            UniformScaleTransform(["x1", "y1"], ["x3", "y3"], 1.5),
        ),
        (
            UniformScaleTransform(["x2", "y2"], ["x3", "y3"], 1.5),
            UniformScaleTransform(["x1", "y1"], ["x2", "y2"], 3.0),
            UniformScaleTransform(["x1", "y1"], ["x3", "y3"], 4.5),
        ),
        (
            UniformScaleTransform(["x2", "y2"], ["x3", "y3"], -0.5),
            ScaleTransform(["x1", "y1"], ["x2", "y2"], [-3.0, 3.0]),
            ScaleTransform(["x1", "y1"], ["x3", "y3"], [1.5, -1.5]),
        ),
        (
            UniformScaleTransform(["x2", "y2"], ["x3", "y3"], 0.5),
            AffineTransform(
                ["x1", "y1"],
                ["x2", "y2"],
                np.array([[1.5, 3.0, 0.0], [-1.5, 3.0, 1.0]], dtype=np.float64),
            ),
            AffineTransform(
                ["x1", "y1"],
                ["x3", "y3"],
                np.array([[0.75, 1.5, 0.0], [-0.75, 1.5, 0.5]], dtype=np.float64),
            ),
        ),
        (
            ScaleTransform(["x2", "y2"], ["x3", "y3"], np.array([1.5, 3.0], dtype=np.float64)),
            IdentityTransform(["x1", "y1"], ["x2", "y2"]),
            ScaleTransform(["x1", "y1"], ["x3", "y3"], np.array([1.5, 3.0], dtype=np.float64)),
        ),
        (
            ScaleTransform(["x2", "y2"], ["x3", "y3"], [1.0, -1.0]),
            UniformScaleTransform(["x1", "y1"], ["x2", "y2"], 1.5),
            ScaleTransform(["x1", "y1"], ["x3", "y3"], [1.5, -1.5]),
        ),
        (
            ScaleTransform(["x2", "y2"], ["x3", "y3"], [1.5, -1.0]),
            ScaleTransform(["x1", "y1"], ["x2", "y2"], [2.0, 1.5]),
            ScaleTransform(["x1", "y1"], ["x3", "y3"], [3.0, -1.5]),
        ),
        (
            ScaleTransform(["x2", "y2"], ["x3", "y3"], [0.5, -0.5]),
            AffineTransform(
                ["x1", "y1"],
                ["x2", "y2"],
                np.array([[1.5, 3.0, 0.0], [-1.5, 3.0, 1.0]], dtype=np.float64),
            ),
            AffineTransform(
                ["x1", "y1"],
                ["x3", "y3"],
                np.array([[0.75, 1.5, 0.0], [0.75, -1.5, -0.5]], dtype=np.float64),
            ),
        ),
        (
            AffineTransform(
                ["x2", "y2"],
                ["x3", "y3"],
                np.array([[1.5, 3.0, 0.0], [-1.5, 3.0, 1.0]], dtype=np.float64),
            ),
            IdentityTransform(["x1", "y1"], ["x2", "y2"]),
            AffineTransform(
                ["x1", "y1"],
                ["x3", "y3"],
                np.array([[1.5, 3.0, 0.0], [-1.5, 3.0, 1.0]], dtype=np.float64),
            ),
        ),
        (
            AffineTransform(
                ["x2", "y2"],
                ["x3", "y3"],
                np.array([[1.5, 3.0, 0.0], [-1.5, 3.0, 1.0]], dtype=np.float64),
            ),
            UniformScaleTransform(["x1", "y1"], ["x2", "y2"], 2.0),
            AffineTransform(
                ["x1", "y1"],
                ["x3", "y3"],
                np.array([[3.0, 6.0, 0.0], [-3.0, 6.0, 1.0]], dtype=np.float64),
            ),
        ),
        (
            AffineTransform(
                ["x2", "y2"],
                ["x3", "y3"],
                np.array([[1.5, 3.0, 0.0], [-1.5, 3.0, 1.0]], dtype=np.float64),
            ),
            ScaleTransform(["x1", "y1"], ["x2", "y2"], [2.0, -2.0]),
            AffineTransform(
                ["x1", "y1"],
                ["x3", "y3"],
                np.array([[3.0, -6.0, 0.0], [-3.0, -6.0, 1.0]], dtype=np.float64),
            ),
        ),
        (
            AffineTransform(
                ["x2", "y2"],
                ["x3", "y3"],
                np.array([[2.0, 0.0, 1.0], [0.0, 4.0, 1.0]], dtype=np.float64),
            ),
            AffineTransform(
                ["x1", "y1"],
                ["x2", "y2"],
                np.array([[1.0, 1.0, -1.0], [0.0, 1.0, 2.0]], dtype=np.float64),
            ),
            AffineTransform(
                ["x1", "y1"],
                ["x3", "y3"],
                np.array([[2.0, 2.0, -1.0], [0.0, 4.0, 9.0]], dtype=np.float64),
            ),
        ),
    ],
    ids=lambda val: type(val).__name__,
)
def test_multiply_tranform(
    transform_a,
    transform_b,
    expected: CoordinateTransform,
):
    result = transform_a @ transform_b
    check_transform_is_equal(result, expected)
