from urllib.parse import urljoin

import numpy as np
import pyarrow as pa
import pytest
import somacore

import tiledbsoma as soma

soma_outgest = pytest.importorskip("tiledbsoma.io.spatial.outgest")
sd = pytest.importorskip("spatialdata")


@pytest.fixture(scope="module")
def sample_2d_data():
    return [
        np.random.randint(0, 255, size=(3, 32, 32), dtype=np.uint8),
        np.random.randint(0, 255, size=(3, 16, 16), dtype=np.uint8),
        np.random.randint(0, 255, size=(3, 8, 8), dtype=np.uint8),
    ]


@pytest.fixture(scope="module")
def sample_multiscale_image_2d(tmp_path_factory, sample_2d_data):
    # Create the multiscale image.
    baseuri = tmp_path_factory.mktemp("export_multiscale_image").as_uri()
    image_uri = urljoin(baseuri, "default")
    with soma.MultiscaleImage.create(
        image_uri,
        type=pa.uint8(),
        coordinate_space=("x_image", "y_image"),
        level_shape=(3, 32, 32),
    ) as image:
        coords = (slice(None), slice(None), slice(None))
        # Create levels.
        l0 = image["level0"]
        l0.write(coords, pa.Tensor.from_numpy(sample_2d_data[0]))

        # Create medium sized downsample.
        l1 = image.add_new_level("level1", shape=(3, 16, 16))
        l1.write(coords, pa.Tensor.from_numpy(sample_2d_data[1]))

        # Create very small downsample and write to it.
        l2 = image.add_new_level("level2", shape=(3, 8, 8))
        l2.write(coords, pa.Tensor.from_numpy(sample_2d_data[2]))
    image2d = soma.MultiscaleImage.open(image_uri)
    return image2d


@pytest.mark.parametrize(
    "level,transform,expected_transformation,expected_transformation_key",
    [
        (
            0,
            None,
            sd.transformations.Identity(),
            "image",
        ),
        (
            2,
            None,
            sd.transformations.Scale([4, 4], ("x", "y")),
            "image",
        ),
        (
            0,
            somacore.IdentityTransform(("x_scene", "y_scene"), ("x_image", "y_image")),
            sd.transformations.Identity(),
            "scene0",
        ),
        (
            2,
            somacore.IdentityTransform(("x_scene", "y_scene"), ("x_image", "y_image")),
            sd.transformations.Scale([4, 4], ("x", "y")),
            "scene0",
        ),
        (
            0,
            somacore.ScaleTransform(
                ("x_scene", "y_scene"), ("x_image", "y_image"), [0.25, 0.5]
            ),
            sd.transformations.Scale([4, 2], ("x", "y")),
            "scene0",
        ),
        (
            2,
            somacore.ScaleTransform(
                ("x_scene", "y_scene"), ("x_image", "y_image"), [0.25, 0.5]
            ),
            sd.transformations.Scale([16, 8], ("x", "y")),
            "scene0",
        ),
        (
            0,
            somacore.AffineTransform(
                ("x_scene", "y_scene"), ("x_image", "y_image"), [[1, 0, 1], [0, 1, 2]]
            ),
            sd.transformations.Affine(
                np.array([[1, 0, -1], [0, 1, -2], [0, 0, 1]]),
                ("x", "y"),
                ("x", "y"),
            ),
            "scene0",
        ),
        (
            2,
            somacore.AffineTransform(
                ("x_scene", "y_scene"), ("x_image", "y_image"), [[1, 0, 1], [0, 1, 2]]
            ),
            sd.transformations.Sequence(
                [
                    sd.transformations.Scale([4, 4], ("x", "y")),
                    sd.transformations.Affine(
                        np.array([[1, 0, -1], [0, 1, -2], [0, 0, 1]]),
                        ("x", "y"),
                        ("x", "y"),
                    ),
                ]
            ),
            "scene0",
        ),
    ],
)
def test_export_image_level_to_spatialdata(
    sample_multiscale_image_2d,
    sample_2d_data,
    level,
    transform,
    expected_transformation,
    expected_transformation_key,
):
    image2d = soma_outgest.to_spatialdata_image(
        sample_multiscale_image_2d,
        level=level,
        key="image",
        scene_id="scene0",
        scene_dim_map={"x_scene": "x", "y_scene": "y"},
        transform=transform,
    )

    assert isinstance(image2d, sd.models.models.DataArray)

    # Validate the model.
    schema = sd.models.get_model(image2d)
    assert schema == sd.models.Image2DModel

    # Check the correct data exists.
    result = image2d.data.compute()
    np.testing.assert_equal(result, sample_2d_data[level])

    # Check the metadata.
    metadata = dict(image2d.attrs)
    assert len(metadata) == 1
    assert metadata["transform"] == {
        expected_transformation_key: expected_transformation
    }


@pytest.mark.parametrize(
    "transform,expected_transformation,expected_transformation_key",
    [
        (
            None,
            [
                sd.transformations.Identity(),
                sd.transformations.Scale([2, 2], ("x", "y")),
                sd.transformations.Scale([4, 4], ("x", "y")),
            ],
            "image",
        ),
        (
            somacore.IdentityTransform(("x_scene", "y_scene"), ("x_image", "y_image")),
            [
                sd.transformations.Identity(),
                sd.transformations.Scale([2, 2], ("x", "y")),
                sd.transformations.Scale([4, 4], ("x", "y")),
            ],
            "scene0",
        ),
        (
            somacore.ScaleTransform(
                ("x_scene", "y_scene"), ("x_image", "y_image"), [0.25, 0.5]
            ),
            [
                sd.transformations.Scale([4, 2], ("x", "y")),
                sd.transformations.Scale([8, 4], ("x", "y")),
                sd.transformations.Scale([16, 8], ("x", "y")),
            ],
            "scene0",
        ),
        (
            somacore.AffineTransform(
                ("x_scene", "y_scene"), ("x_image", "y_image"), [[1, 0, 1], [0, 1, 2]]
            ),
            [
                sd.transformations.Affine(
                    np.array([[1, 0, -1], [0, 1, -2], [0, 0, 1]]),
                    ("x", "y"),
                    ("x", "y"),
                ),
                sd.transformations.Sequence(
                    [
                        sd.transformations.Scale([2, 2], ("x", "y")),
                        sd.transformations.Affine(
                            np.array([[1, 0, -1], [0, 1, -2], [0, 0, 1]]),
                            ("x", "y"),
                            ("x", "y"),
                        ),
                    ]
                ),
                sd.transformations.Sequence(
                    [
                        sd.transformations.Scale([4, 4], ("x", "y")),
                        sd.transformations.Affine(
                            np.array([[1, 0, -1], [0, 1, -2], [0, 0, 1]]),
                            ("x", "y"),
                            ("x", "y"),
                        ),
                    ]
                ),
            ],
            "scene0",
        ),
    ],
)
def test_export_full_image_to_spatialdata(
    sample_multiscale_image_2d,
    sample_2d_data,
    transform,
    expected_transformation,
    expected_transformation_key,
):
    image2d = soma_outgest.to_spatialdata_multiscale_image(
        sample_multiscale_image_2d,
        key="image",
        scene_id="scene0",
        scene_dim_map={"x_scene": "x", "y_scene": "y"},
        transform=transform,
    )

    assert isinstance(image2d, sd.models.models.DataTree)

    # Validate the model.
    schema = sd.models.get_model(image2d)
    assert schema == sd.models.Image2DModel

    # Check the correct data exists.
    for index in range(3):
        data_array = image2d[f"scale{index}"]["image"]

        # Check data.
        result = data_array.data.compute()
        np.testing.assert_equal(result, sample_2d_data[index])

        # Check the metadata.
        metadata = dict(data_array.attrs)
        assert len(metadata) == 1
        assert metadata["transform"] == {
            expected_transformation_key: expected_transformation[index]
        }
