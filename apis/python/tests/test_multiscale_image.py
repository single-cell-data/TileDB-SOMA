import functools
from urllib.parse import urljoin

import numpy as np
import pyarrow as pa
import pytest

import tiledbsoma as soma
from tiledbsoma import IdentityTransform, ScaleTransform


def test_multiscale_image_bad_create(tmp_path):
    baseuri = urljoin(f"{tmp_path.as_uri()}/", "bad_create")

    # Invalid datetype.
    with pytest.raises(ValueError):
        soma.MultiscaleImage.create(
            baseuri,
            type=pa.string(),
            coordinate_space=("x", "y", "y"),
            level_shape=(3, 64, 64),
        )

    # Repeated axis names.
    with pytest.raises(ValueError):
        soma.MultiscaleImage.create(
            baseuri,
            type=pa.uint8(),
            coordinate_space=("x", "y", "y"),
            level_shape=(3, 64, 64),
        )

    # Repeated data axis order, missing "soma_channel".
    with pytest.raises(ValueError):
        soma.MultiscaleImage.create(
            baseuri,
            type=pa.uint8(),
            data_axis_order=("x", "x", "y"),
            level_shape=(3, 64, 64),
        )

    # Repeated data axis order, bad length.
    with pytest.raises(ValueError):
        soma.MultiscaleImage.create(
            baseuri,
            type=pa.uint8(),
            data_axis_order=("soma_channel", "x", "x", "y"),
            level_shape=(3, 64, 64),
        )

    # Invalid data axis order.
    with pytest.raises(ValueError):
        soma.MultiscaleImage.create(
            baseuri,
            type=pa.uint8(),
            data_axis_order=("soma_channel", "y", "bad-axis-name"),
            level_shape=(3, 64, 64),
        )

    # Data order has `soma_channel` when not expected
    with pytest.raises(ValueError):
        soma.MultiscaleImage.create(
            baseuri,
            type=pa.uint8(),
            data_axis_order=("soma_channel", "y", "x"),
            level_shape=(3, 64, 64),
            has_channel_axis=False,
        )


def test_multiscale_basic(tmp_path):
    baseuri = urljoin(f"{tmp_path.as_uri()}/", "basic_read")
    image_uri = urljoin(baseuri, "default")

    # Create the multiscale image.
    with soma.MultiscaleImage.create(
        image_uri,
        type=pa.uint8(),
        level_shape=(128, 64),
        has_channel_axis=False,
    ) as image:

        # Add medium sized downsample.
        image.add_new_level("level1", shape=(64, 32))

        # Add very small downsample and write to it.
        level2 = image.add_new_level("level2", shape=(8, 4))
        level2_data = pa.Tensor.from_numpy(np.arange(32, dtype=np.uint8).reshape(8, 4))
        level2.write((slice(None), slice(None)), level2_data)

    # Open for reading and check metadata.
    with soma.MultiscaleImage.open(image_uri, mode="r") as image:

        # Check the base properties for the image.
        assert image.level_shape(0) == (128, 64)
        assert image.nchannels == 1
        assert image.data_axis_order == ("y", "x")

        # Check coordinate space.
        coord_space = image.coordinate_space
        assert len(coord_space) == 2
        assert coord_space.axis_names == ("x", "y")

        # Check the number of levels and level properties.
        expected_shapes = [(128, 64), (64, 32), (8, 4)]
        assert image.level_count == 3
        for index, shape in enumerate(expected_shapes):
            assert shape == expected_shapes[index]
            assert image.level_shape(index) == expected_shapes[index]

        # Check the spatial version encoding was written.
        assert (
            image.metadata[soma._constants.SOMA_SPATIAL_VERSION_METADATA_KEY]
            == soma._constants.SOMA_SPATIAL_ENCODING_VERSION
        )

        # Check the levels mapping.
        levels = image.levels()
        assert len(levels) == 3
        for key, (uri, shape) in levels.items():
            assert soma.DenseNDArray.exists(uri)
            level_image = soma.DenseNDArray.open(uri)
            assert level_image.shape == image.level_shape(key)
            assert shape == image.level_shape(key)

        # Check the level_uri function.
        for index in range(3):
            uri = image.level_uri(index)
            print(uri)
            assert soma.DenseNDArray.exists(uri)
            level_image = soma.DenseNDArray.open(uri)
            assert level_image.shape == expected_shapes[index]

        for index, key in enumerate(["level0", "level1", "level2"]):
            uri = image.level_uri(key)
            assert soma.DenseNDArray.exists(uri)
            level_image = soma.DenseNDArray.open(uri)
            assert level_image.shape == expected_shapes[index]

        # Check a basic read
        assert level2_data == image.read_spatial_region(2).data

        # Check transform to and from levels
        to_level = image.get_transform_to_level
        from_level = image.get_transform_from_level

        assert isinstance(to_level(0), IdentityTransform)
        assert isinstance(from_level(0), IdentityTransform)
        assert np.array_equal(to_level(0).scale_factors, [1, 1])
        assert np.array_equal(to_level(1).scale_factors, [0.5, 0.5])
        assert np.array_equal(to_level(2).scale_factors, [0.0625, 0.0625])
        assert np.array_equal(to_level("level0").scale_factors, [1, 1])
        assert np.array_equal(to_level("level1").scale_factors, [0.5, 0.5])
        assert np.array_equal(to_level("level2").scale_factors, [0.0625, 0.0625])

        # oob
        with pytest.raises(IndexError):
            to_level(3).scale_factors

        # dne
        with pytest.raises(KeyError):
            to_level("level3").scale_factors

        assert np.array_equal(from_level(0).scale_factors, [1, 1])
        assert np.array_equal(from_level(1).scale_factors, [2, 2])
        assert np.array_equal(from_level(2).scale_factors, [16, 16])

    # Ensure it cannot be opened by another type
    with pytest.raises(soma.SOMAError):
        soma.DataFrame.open(image_uri)

    with pytest.raises(soma.SOMAError):
        soma.SparseNDArray.open(image_uri)

    with pytest.raises(soma.SOMAError):
        soma.DenseNDArray.open(image_uri)

    with pytest.raises(soma.SOMAError):
        soma.PointCloudDataFrame.open(image_uri)

    with pytest.raises(soma.SOMAError):
        soma.Collection.open(image_uri)

    with pytest.raises(soma.SOMAError):
        soma.Experiment.open(image_uri)

    with pytest.raises(soma.SOMAError):
        soma.Measurement.open(image_uri)

    with pytest.raises(soma.SOMAError):
        soma.Scene.open(image_uri)


class TestSimpleMultiscale2D:

    @pytest.fixture(scope="class")
    def image_uri(self, tmp_path_factory):
        """Create a multiscale image and return the path."""
        # Create the multiscale image.
        baseuri = tmp_path_factory.mktemp("multiscale_image").as_uri()
        image_uri = urljoin(baseuri, "simple2d")
        with soma.MultiscaleImage.create(
            image_uri,
            type=pa.uint8(),
            level_shape=(1, 9, 8),
        ) as image:
            coords = (slice(None), slice(None), slice(None))
            # Create levels.
            l0 = image["level0"]
            l0.write(
                coords,
                pa.Tensor.from_numpy(np.arange(72, dtype=np.uint8).reshape(1, 9, 8)),
            )

            # Create medium sized downsample.
            l1 = image.add_new_level("level1", shape=(1, 6, 4))
            l1.write(
                coords,
                pa.Tensor.from_numpy(
                    10 * np.arange(24, dtype=np.uint8).reshape(1, 6, 4)
                ),
            )

            # Create very small downsample and write to it.
            l2 = image.add_new_level("level2", shape=(1, 3, 2))
            l2.write(
                coords,
                pa.Tensor.from_numpy(
                    100 * np.arange(6, dtype=np.uint8).reshape(1, 3, 2)
                ),
            )
        return image_uri

    @pytest.mark.parametrize(
        ("level", "region", "kwargs", "expected_data", "expected_transform"),
        [
            pytest.param(
                2,
                None,
                {},
                100 * np.arange(6, dtype=np.uint8).reshape(1, 3, 2),
                ScaleTransform(("x", "y"), ("x", "y"), [4, 3]),
                id="Level 2, full region, no transform",
            ),
        ],
    )
    def test_read_spatial_region(
        self, image_uri, level, region, kwargs, expected_data, expected_transform
    ):
        with soma.MultiscaleImage.open(image_uri) as image:
            result = image.read_spatial_region(level=level, region=region, **kwargs)
        actual_data = result.data.to_numpy()

        # Check data
        np.testing.assert_array_equal(actual_data, expected_data)

        # Check transform
        actual_transform = result.coordinate_transform
        assert actual_transform.input_axes == expected_transform.input_axes
        assert actual_transform.output_axes == expected_transform.output_axes
        assert isinstance(actual_transform, type(expected_transform))
        np.testing.assert_array_almost_equal(
            actual_transform.augmented_matrix,
            expected_transform.augmented_matrix,
            decimal=8,
        )


def create_multiscale(baseuri, coord_space, axis_order, has_channel_axis, shapes):
    image_uri = urljoin(baseuri, "default")
    with soma.MultiscaleImage.create(
        image_uri,
        type=pa.uint8(),
        coordinate_space=coord_space,
        data_axis_order=axis_order,
        has_channel_axis=has_channel_axis,
        level_shape=shapes[0],
    ) as image:
        for i in range(1, len(shapes)):
            image.add_new_level(f"level{i}", shape=shapes[i])
    return image_uri


@pytest.mark.parametrize(
    "coord_space, axis_order, has_channel_axis, shapes, expected_scale_factors",
    [
        [
            ("x", "y", "z"),
            ("soma_channel", "z", "y", "x"),
            True,
            ((128, 64, 32, 16), (128, 32, 16, 8), (128, 16, 8, 4), (128, 4, 2, 1)),
            ([1, 1, 1], [2, 2, 2], [4, 4, 4], [16, 16, 16]),
        ],
        [
            ("x", "y", "z"),
            ("soma_channel", "z", "y", "x"),
            True,
            ((128, 64, 32, 16), (128, 32, 16, 8), (128, 16, 8, 4)),
            ([1, 1, 1], [2, 2, 2], [4, 4, 4]),
        ],
        [
            ("x", "y", "z"),
            ("soma_channel", "z", "y", "x"),
            True,
            ((64, 64, 64, 64), (64, 64, 64, 64)),
            ([1, 1, 1], [1, 1, 1]),
        ],
        [
            ("x", "y", "z"),
            ("soma_channel", "z", "y", "x"),
            True,
            ((64, 32, 16, 8), (64, 16, 8, 4)),
            ([1, 1, 1], [2, 2, 2]),
        ],
        [
            ("x", "y", "z"),
            ("soma_channel", "z", "y", "x"),
            True,
            ((128, 64, 32, 16), (128, 32, 32, 8), (128, 16, 16, 4)),
            ([1, 1, 1], [2, 1, 2], [4, 2, 4]),
        ],
        [
            ("x", "y"),
            ("soma_channel", "y", "x"),
            True,
            ((128, 64, 32), (128, 32, 16), (128, 16, 8)),
            ([1, 1], [2, 2], [4, 4]),
        ],
        [
            ("x", "y"),
            ("soma_channel", "y", "x"),
            True,
            ((128, 128, 128), (128, 128, 128)),
            ([1, 1], [1, 1]),
        ],
        [
            ("x", "y"),
            ("y", "x"),
            False,
            ((128, 128), (128, 128)),
            ([1, 1], [1, 1]),
        ],
        [
            ("x", "y"),
            ("y", "x"),
            False,
            ((128, 64), (64, 32)),
            ([1, 1], [2, 2]),
        ],
        [
            ("x", "y"),
            ("y", "x"),
            False,
            ((60, 30), (30, 6)),
            ([1, 1], [5, 2]),
        ],
    ],
)
def test_multiscale_with_axis_names(
    tmp_path, coord_space, axis_order, has_channel_axis, shapes, expected_scale_factors
):
    baseuri = urljoin(f"{tmp_path.as_uri()}/", "test_multiscale_with_axis_names")
    image_uri = create_multiscale(
        baseuri, coord_space, axis_order, has_channel_axis, shapes
    )

    with soma.MultiscaleImage.open(image_uri, mode="r") as image:
        assert image.level_count == len(shapes)

        # TODO Check channels

        for index, shape in enumerate(shapes):
            assert shape == image.level_shape(index)

            # Check transform to levels.
            expected = 1 / np.array(expected_scale_factors[index])
            actual = image.get_transform_to_level(index).scale_factors
            assert np.array_equal(actual, expected)
            assert np.array_equal(actual, expected)

            # Check transform from levels
            expected = expected_scale_factors[index]
            actual = image.get_transform_from_level(index).scale_factors
            assert np.array_equal(actual, expected)
            actual = image.get_transform_from_level(f"level{index}").scale_factors
            assert np.array_equal(actual, expected)


@pytest.mark.parametrize(
    "shapes, data_axis_order, region, expected_coords, scale_factors",
    [
        # full region
        (
            ((3, 64, 32), (3, 32, 16), (3, 16, 8)),
            ("soma_channel", "y", "x"),
            None,
            None,
            ([1, 1], [2, 2], [4, 4]),
        ),
        (
            ((64, 32, 3), (32, 16, 3), (16, 8, 3)),
            ("y", "x", "soma_channel"),
            None,
            None,
            ([1, 1], [2, 2], [4, 4]),
        ),
        (
            ((64, 3, 32), (32, 3, 16), (16, 3, 8)),
            ("x", "soma_channel", "y"),
            None,
            None,
            ([1, 1], [2, 2], [4, 4]),
        ),
        (
            ((3, 128, 128), (3, 128, 128)),
            ("soma_channel", "y", "x"),
            None,
            None,
            ([1, 1], [1, 1]),
        ),
        (
            ((3, 128, 64), (3, 64, 32)),
            ("soma_channel", "y", "x"),
            None,
            None,
            ([1, 1], [2, 2]),
        ),
        (
            ((3, 60, 30), (3, 30, 6)),
            ("soma_channel", "y", "x"),
            None,
            None,
            ([1, 1], [2, 5]),
        ),
        (
            ((3, 1, 1),),
            ("soma_channel", "y", "x"),
            None,
            None,
            ([1, 1]),
        ),
        # partial region
        pytest.param(
            ((3, 128, 64), (3, 64, 32)),
            ("soma_channel", "y", "x"),
            (0, 0, 20, 30),
            (
                (slice(None), slice(0, 30), slice(0, 20)),
                (slice(None), slice(0, 15), slice(0, 10)),
            ),
            ([1, 1], [2, 2]),
            id="partial region - channel first",
        ),
        pytest.param(
            ((64, 3, 128), (32, 3, 64)),
            ("x", "soma_channel", "y"),
            (0, 0, 20, 30),
            (
                (slice(0, 20), slice(None), slice(0, 30)),
                (slice(0, 10), slice(None), slice(0, 15)),
            ),
            ([1, 1], [2, 2]),
            id="partial region - channel middle",
        ),
        pytest.param(
            ((128, 64, 3), (64, 32, 3)),
            ("y", "x", "soma_channel"),
            (0, 0, 20, 30),
            (
                (slice(0, 30), slice(0, 20), slice(None)),
                (slice(0, 15), slice(0, 10), slice(None)),
            ),
            ([1, 1], [2, 2]),
            id="partial region - channel last",
        ),
        pytest.param(
            ((3, 64, 32), (3, 32, 16)),
            ("soma_channel", "y", "x"),
            (0, 0, 16, 10),
            (
                (slice(None), slice(0, 10), slice(0, 16)),
                (slice(None), slice(0, 5), slice(0, 8)),
            ),
            ([1, 1], [2, 2]),
            id="partial region  - small",
        ),
    ],
)
def test_multiscale_2d_read_region_with_channel(
    tmp_path, shapes, data_axis_order, region, expected_coords, scale_factors
):
    baseuri = urljoin(f"{tmp_path.as_uri()}/", "test_multiscale_read_region")
    image_uri = create_multiscale(baseuri, ("x", "y"), data_axis_order, True, shapes)

    with soma.MultiscaleImage.open(image_uri, mode="w") as image:
        for i, shape in enumerate(shapes):
            data = np.arange(shape[0] * shape[1] * shape[2], dtype=np.uint8).reshape(
                shape
            )
            image[f"level{i}"].write(
                (slice(None), slice(None)), pa.Tensor.from_numpy(data)
            )

    with soma.MultiscaleImage.open(image_uri, mode="r") as image:
        for i, shape in enumerate(shapes):
            actual_data = image.read_spatial_region(i, region=region).data
            if expected_coords is None:
                expected_data = image[f"level{i}"].read()
            else:
                expected_data = image[f"level{i}"].read(expected_coords[i])
            assert np.array_equal(actual_data, expected_data)


@pytest.mark.parametrize(
    "shapes, region, scale_factors",
    [
        # full region
        (
            ((64, 32), (32, 16), (16, 8)),
            None,
            ([1, 1], [2, 2], [4, 4]),
        ),
        (
            ((128, 128), (128, 128)),
            None,
            ([1, 1], [1, 1]),
        ),
        (
            ((128, 64), (64, 32)),
            None,
            ([1, 1], [2, 2]),
        ),
        (
            ((60, 30), (30, 6)),
            None,
            ([1, 1], [2, 5]),
        ),
        (
            ((1, 1),),
            None,
            ([1, 1]),
        ),
        # partial region
        (
            ((128, 64), (64, 32)),
            (0, 0, 20, 30),
            ([1, 1], [2, 2]),
        ),
        (
            ((64, 32), (32, 16)),
            (0, 0, 16, 10),
            ([1, 1], [2, 2]),
        ),
    ],
)
def test_multiscale_2d_read_region_no_channel(tmp_path, shapes, region, scale_factors):
    baseuri = urljoin(f"{tmp_path.as_uri()}/", "test_multiscale_read_region")
    image_uri = create_multiscale(baseuri, ("x", "y"), ("y", "x"), False, shapes)

    with soma.MultiscaleImage.open(image_uri, mode="w") as image:
        for i, shape in enumerate(shapes):
            data = np.arange(shape[0] * shape[1], dtype=np.uint8).reshape(shape)
            image[f"level{i}"].write(
                (slice(None), slice(None)), pa.Tensor.from_numpy(data)
            )

    with soma.MultiscaleImage.open(image_uri, mode="r") as image:
        for i, shape in enumerate(shapes):
            actual_data = image.read_spatial_region(i, region=region).data
            if region is None:
                expected_data = image[f"level{i}"].read()
            else:
                expected_data = image[f"level{i}"].read(
                    coords=(
                        slice(
                            region[1],
                            region[3] // scale_factors[i][0],
                        ),
                        slice(
                            region[0],
                            region[2] // scale_factors[i][1],
                        ),
                    )
                )
            assert np.array_equal(actual_data, expected_data)


@pytest.mark.skip("reading 3D regions not supported yet")
@pytest.mark.parametrize(
    "shapes, region, scale_factors",
    [
        # full region
        (
            ((64, 32, 16), (32, 16, 8), (16, 8, 4), (4, 2, 1)),
            None,
            ([1, 1, 1], [2, 2, 2], [4, 4, 4], [16, 16, 16]),
        ),
        (
            ((64, 32, 16), (32, 16, 8), (16, 8, 4)),
            None,
            ([1, 1, 1], [2, 2, 2], [4, 4, 4]),
        ),
        (
            ((64, 64, 64), (64, 64, 64)),
            None,
            ([1, 1, 1], [1, 1, 1]),
        ),
        (
            ((32, 16, 8), (16, 8, 4)),
            None,
            ([1, 1, 1], [2, 2, 2]),
        ),
        (
            ((64, 32, 16), (32, 32, 8), (16, 16, 4)),
            None,
            ([1, 1, 1], [2, 1, 2], [4, 2, 2]),
        ),
        # partial region
        (
            ((64, 32, 16), (32, 16, 8)),
            (0, 0, 0, 16, 16, 8),
            ([1, 1, 1], [2, 2, 2]),
        ),
        (
            ((64, 64, 64), (64, 64, 64)),
            (0, 0, 0, 48, 48, 48),
            ([1, 1, 1], [1, 1, 1]),
        ),
    ],
)
def test_multiscale_3d_read_region_no_channel(tmp_path, shapes, region, scale_factors):
    baseuri = urljoin(f"{tmp_path.as_uri()}/", "test_multiscale_read_region")
    image_uri = create_multiscale(
        baseuri, ("x", "y", "z"), ("z", "y", "x"), False, shapes
    )

    with soma.Collection.open(image_uri, mode="w") as image:
        for i, shape in enumerate(shapes):
            size = functools.reduce(lambda x, y: x * y, shape)
            data = np.arange(size, dtype=np.uint8).reshape(*shape)
            image[f"level{i}"].write(
                (slice(None), slice(None), slice(None)), pa.Tensor.from_numpy(data)
            )

    with soma.MultiscaleImage.open(image_uri, mode="r") as image:
        for i, shape in enumerate(shapes):
            if region is None:
                actual_data = image.read_spatial_region(i).data
                expected_data = np.arange(
                    functools.reduce(lambda x, y: x * y, shape), dtype=np.uint8
                ).reshape(*shape)
            else:
                actual_data = image.read_spatial_region(i, region=region).data
                expected_data = image[f"level{i}"].read(
                    coords=(
                        slice(region[2], region[5] // scale_factors[i][2]),
                        slice(region[1], region[4] // scale_factors[i][1]),
                        slice(region[0], region[3] // scale_factors[i][0]),
                    )
                )
            assert np.array_equal(actual_data, expected_data)
