import functools
from urllib.parse import urljoin

import numpy as np
import pyarrow as pa
import pytest

import tiledbsoma as soma
from tiledbsoma import ScaleTransform


def test_multiscale_image_bad_create(tmp_path):
    baseuri = urljoin(f"{tmp_path.as_uri()}/", "bad_create")

    # Invalid datetype.
    with pytest.raises(ValueError):
        soma.MultiscaleImage.create(
            baseuri,
            type=pa.string(),
            axis_names=("x", "y", "y"),
            reference_level_shape=(3, 64, 64),
        )

    # Repeated axis names.
    with pytest.raises(ValueError):
        soma.MultiscaleImage.create(
            baseuri,
            type=pa.uint8(),
            axis_names=("x", "y", "y"),
            reference_level_shape=(3, 64, 64),
        )

    # Repeated axis types.
    with pytest.raises(ValueError):
        soma.MultiscaleImage.create(
            baseuri,
            type=pa.uint8(),
            axis_types=("height", "height", "width"),
            reference_level_shape=(3, 64, 64),
        )

    # Invalid axis type.
    with pytest.raises(ValueError):
        soma.MultiscaleImage.create(
            baseuri,
            type=pa.uint8(),
            axis_types=("c", "height", "width"),
            reference_level_shape=(3, 64, 64),
        )

    # Mismatch in number of axis names and reference shape.
    with pytest.raises(ValueError):
        soma.MultiscaleImage.create(
            baseuri,
            type=pa.uint8(),
            axis_types=("x", "y", "y"),
            reference_level_shape=(3, 64, 64),
        )

    # Mismatch in number of axis names and axis types.
    with pytest.raises(ValueError):
        soma.MultiscaleImage.create(
            baseuri,
            type=pa.uint8(),
            reference_level_shape=(64, 64),
            axis_names=("y", "x"),
            axis_types=("channels",),
        )


def test_multiscale_basic(tmp_path):
    baseuri = urljoin(f"{tmp_path.as_uri()}/", "basic_read")
    image_uri = urljoin(baseuri, "default")

    # Create the multiscale image.
    with soma.MultiscaleImage.create(
        image_uri,
        type=pa.uint8(),
        reference_level_shape=(128, 64),
        axis_names=("y", "x"),
        axis_types=("height", "width"),
    ) as image:

        # Create reference level.
        image.add_new_level("level0", shape=(128, 64))

        # Create medium sized downsample.
        image.add_new_level("level1", shape=(64, 32))

        # Create very small downsample and write to it.
        level2 = image.add_new_level("level2", shape=(8, 4))
        level2_data = pa.Tensor.from_numpy(np.arange(32, dtype=np.uint8).reshape(8, 4))
        level2.write((slice(None), slice(None)), level2_data)

    # Open for reading and check metadata.
    with soma.MultiscaleImage.open(image_uri, mode="r") as image:

        # Check the base properties for the image.
        assert image.axis_names == ("y", "x")
        base_props = image.reference_level_properties
        assert base_props.name == "reference_level"
        assert base_props.shape == (128, 64)
        assert base_props.height == 128
        assert base_props.width == 64
        assert base_props.nchannels is None
        assert base_props.depth is None
        assert base_props.image_type == "YX"

        # Check coordinate space.
        coord_space = image.coordinate_space
        assert len(coord_space) == 2
        assert coord_space.axis_names == ("x", "y")

        # Check the number of levels and level properties.
        assert image.level_count == 3
        for index, shape in enumerate([(128, 64), (64, 32), (8, 4)]):
            props = image.level_properties(index)
            assert props == image.level_properties(props.name)
            assert props.nchannels is None
            assert props.depth is None
            assert props.image_type == "YX"
            assert props.name == f"level{index}"
            assert props.shape == shape
            assert props.height == shape[0]
            assert props.width == shape[1]

        # Check a basic read
        assert level2_data == image.read_spatial_region(2).data

        # Check transform to and from levels
        to_level = image.get_transform_to_level
        from_level = image.get_transform_from_level

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


class TestSimpleMultiscale2D:

    @pytest.fixture(scope="class")
    def image_uri(self, tmp_path_factory):
        """Create a multiscale image and return the path."""
        # Create the multiscale image.
        baseuri = tmp_path_factory.mktemp("multiscale_image").as_uri()
        image_uri = urljoin(baseuri, "default")
        with soma.MultiscaleImage.create(
            image_uri,
            type=pa.uint8(),
            reference_level_shape=(1, 9, 8),
            axis_names=("c", "y", "x"),
            axis_types=("channel", "height", "width"),
        ) as image:
            coords = (slice(None), slice(None), slice(None))
            # Create levels.
            l0 = image.add_new_level("level0", shape=(1, 9, 8))
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


def create_multiscale(baseuri, axis_names, axis_types, shapes):
    image_uri = urljoin(baseuri, "default")
    with soma.MultiscaleImage.create(
        image_uri,
        type=pa.uint8(),
        axis_names=axis_names,
        axis_types=axis_types,
        reference_level_shape=shapes[0],
    ) as image:
        for i in range(len(shapes)):
            image.add_new_level(f"level{i}", shape=shapes[i])
    return image_uri


@pytest.mark.parametrize(
    "axis_names, axis_types, shapes, expected_scale_factors",
    [
        [
            ("C", "Z", "Y", "X"),
            ("channel", "depth", "height", "width"),
            ((128, 64, 32, 16), (128, 32, 16, 8), (128, 16, 8, 4), (128, 4, 2, 1)),
            ([1, 1, 1], [2, 2, 2], [4, 4, 4], [16, 16, 16]),
        ],
        [
            ("C", "Z", "Y", "X"),
            ("channel", "depth", "height", "width"),
            ((128, 64, 32, 16), (128, 32, 16, 8), (128, 16, 8, 4)),
            ([1, 1, 1], [2, 2, 2], [4, 4, 4]),
        ],
        [
            ("C", "Z", "Y", "X"),
            ("channel", "depth", "height", "width"),
            ((64, 64, 64, 64), (64, 64, 64, 64)),
            ([1, 1, 1], [1, 1, 1]),
        ],
        [
            ("C", "Z", "Y", "X"),
            ("channel", "depth", "height", "width"),
            ((64, 32, 16, 8), (64, 16, 8, 4)),
            ([1, 1, 1], [2, 2, 2]),
        ],
        [
            ("C", "Z", "Y", "X"),
            ("channel", "depth", "height", "width"),
            ((128, 64, 32, 16), (128, 32, 32, 8), (128, 16, 16, 4)),
            ([1, 1, 1], [2, 1, 2], [4, 2, 4]),
        ],
        [
            ("C", "Y", "X"),
            ("channel", "height", "width"),
            ((128, 64, 32), (128, 32, 16), (128, 16, 8)),
            ([1, 1], [2, 2], [4, 4]),
        ],
        [
            ("C", "Y", "X"),
            ("channel", "height", "width"),
            ((128, 128, 128), (128, 128, 128)),
            ([1, 1], [1, 1]),
        ],
        [
            ("Y", "X"),
            ("height", "width"),
            ((128, 128), (128, 128)),
            ([1, 1], [1, 1]),
        ],
        [
            ("Y", "X"),
            ("height", "width"),
            ((128, 64), (64, 32)),
            ([1, 1], [2, 2]),
        ],
        [
            ("Y", "X"),
            ("height", "width"),
            ((60, 30), (30, 6)),
            ([1, 1], [5, 2]),
        ],
    ],
)
def test_multiscale_with_axis_names(
    tmp_path, axis_names, axis_types, shapes, expected_scale_factors
):
    baseuri = urljoin(f"{tmp_path.as_uri()}/", "test_multiscale_with_axis_names")
    image_uri = create_multiscale(baseuri, axis_names, axis_types, shapes)

    with soma.MultiscaleImage.open(image_uri, mode="r") as image:
        base_props = image.reference_level_properties
        assert base_props.shape == shapes[0]
        assert image.level_count == len(shapes)

        for index, shape in enumerate(shapes):
            props = image.level_properties(index)
            assert props.name == f"level{index}"
            assert props == image.level_properties(props.name)
            assert props.image_type == "".join(axis_names)
            assert props.shape == shape

            for i, axis_type in enumerate(axis_types):
                if axis_type == "channel":
                    assert getattr(props, "nchannels") == shape[i]
                else:
                    assert getattr(props, axis_type) == shape[i]

            # Check transform to and from levels
            assert np.array_equal(
                image.get_transform_to_level(index).scale_factors,
                1 / np.array(expected_scale_factors[index]),
            )
            assert np.array_equal(
                image.get_transform_to_level(f"level{index}").scale_factors,
                1 / np.array(expected_scale_factors[index]),
            )
            assert np.array_equal(
                image.get_transform_from_level(index).scale_factors,
                expected_scale_factors[index],
            )
            assert np.array_equal(
                image.get_transform_from_level(f"level{index}").scale_factors,
                expected_scale_factors[index],
            )


@pytest.mark.parametrize(
    "shapes, region",
    [
        # full region
        (
            ((64, 32), (32, 16), (16, 8)),
            None,
        ),
        (
            ((128, 128), (128, 128)),
            None,
        ),
        (
            ((128, 64), (64, 32)),
            None,
        ),
        (
            ((60, 30), (30, 6)),
            None,
        ),
        (
            ((1, 1),),
            None,
        ),
        # partial subregion
        (
            ((128, 64), (64, 32)),
            (0, 0, 20, 30),
        ),
        (
            ((64, 32), (32, 16)),
            (0, 0, 16, 10),
        ),
    ],
)
def test_multiscale_2d_read_region(tmp_path, shapes, region):
    baseuri = urljoin(f"{tmp_path.as_uri()}/", "test_multiscale_read_region")
    image_uri = create_multiscale(baseuri, ("Y", "X"), ("height", "width"), shapes)

    with soma.Collection.open(image_uri, mode="w") as image:
        for i, shape in enumerate(shapes):
            data = np.arange(shape[0] * shape[1], dtype=np.uint8).reshape(shape)
            image[f"level{i}"].write(
                (slice(None), slice(None)), pa.Tensor.from_numpy(data)
            )

    with soma.MultiscaleImage.open(image_uri, mode="r") as image:
        for i, shape in enumerate(shapes):
            actual_data = image.read_spatial_region(i, region=region).data
            expected_data = np.arange(shape[0] * shape[1], dtype=np.uint8).reshape(
                shape
            )
            if region is not None:
                expected_data = expected_data[
                    region[1] : region[3] + 1, region[0] : region[2] + 1
                ]
            print("expected_data:", expected_data.size)
            print("actual_data:", actual_data.to_numpy().size)
            # assert (expected_data == actual_data).all()


@pytest.mark.skip("reading 3D regions not supported yet")
@pytest.mark.parametrize(
    "shapes, region",
    [
        (
            ((64, 32, 16), (32, 16, 8), (16, 8, 4), (4, 2, 1)),
            (slice(None), slice(None)),
        ),
        (
            ((64, 32, 16), (32, 16, 8), (16, 8, 4)),
            (slice(None), slice(None)),
        ),
        (
            ((64, 64, 64), (64, 64, 64)),
            (slice(None), slice(None)),
        ),
        (
            ((32, 16, 8), (16, 8, 4)),
            (slice(None), slice(None)),
        ),
        (
            ((64, 32, 16), (32, 32, 8), (16, 16, 4)),
            (slice(None), slice(None)),
        ),
    ],
)
def test_multiscale_3d_read_region(tmp_path, shapes, region):
    baseuri = urljoin(f"{tmp_path.as_uri()}/", "test_multiscale_read_region")
    image_uri = create_multiscale(
        baseuri, ("Z", "Y", "X"), ("depth", "height", "width"), shapes
    )

    with soma.Collection.open(image_uri, mode="w") as image:
        for i, shape in enumerate(shapes):
            size = functools.reduce(lambda x, y: x * y, shape)
            data = np.arange(size, dtype=np.uint8).reshape(*shape)
            image[f"level{i}"].write(region, pa.Tensor.from_numpy(data))

    with soma.MultiscaleImage.open(image_uri, mode="r") as image:
        for i, shape in enumerate(shapes):
            size = functools.reduce(lambda x, y: x * y, shape)
            expected_data = np.arange(size, dtype=np.uint8).reshape(*shape)
            assert np.array_equal(image.read_spatial_region(i).data, expected_data)
