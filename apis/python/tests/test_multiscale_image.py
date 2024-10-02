from urllib.parse import urljoin

import numpy as np
import pyarrow as pa
import pytest

import tiledbsoma as soma


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


def test_multiscale_basic_no_channels(tmp_path):
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
        data = pa.Tensor.from_numpy(np.arange(32, dtype=np.uint8))
        level2.write((slice(None), slice(None)), data)

    # Open for reading and check metadata.
    with soma.MultiscaleImage.open(image_uri, mode="r") as image:

        # Check the base properties for the image.
        # assert image.axis_names == ("y", "x")
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
            assert props.nchannels is None
            assert props.depth is None
            assert props.image_type == "YX"
            assert props.name == f"level{index}"
            assert props.shape == shape
            assert props.height == shape[0]
            assert props.width == shape[1]
