from urllib.parse import urljoin

import numpy as np
import pyarrow as pa
import pytest
import somacore

import tiledbsoma as soma

soma_outgest = pytest.importorskip("tiledbsoma.experimental.outgest")
sd = pytest.importorskip("spatialdata")
xr = pytest.importorskip("xarray")


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
    "level,transform,expected_transformation",
    [
        (
            0,
            somacore.IdentityTransform(("x_scene", "y_scene"), ("x_image", "y_image")),
            sd.transformations.Identity(),
        ),
        (
            2,
            somacore.IdentityTransform(("x_scene", "y_scene"), ("x_image", "y_image")),
            sd.transformations.Scale([4, 4], ("x", "y")),
        ),
        (
            0,
            somacore.ScaleTransform(
                ("x_scene", "y_scene"), ("x_image", "y_image"), [0.25, 0.5]
            ),
            sd.transformations.Scale([4, 2], ("x", "y")),
        ),
        (
            2,
            somacore.ScaleTransform(
                ("x_scene", "y_scene"), ("x_image", "y_image"), [0.25, 0.5]
            ),
            sd.transformations.Scale([16, 8], ("x", "y")),
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
        ),
    ],
)
def test_export_image_level_to_spatial_data(
    sample_multiscale_image_2d,
    sample_2d_data,
    level,
    transform,
    expected_transformation,
):
    image2d = soma_outgest.to_spatial_data_image(
        sample_multiscale_image_2d,
        level=level,
        scene_id="scene0",
        scene_dim_map={"x_scene": "x", "y_scene": "y"},
        transform=transform,
    )

    assert isinstance(image2d, xr.DataArray)

    # Validate the model.
    schema = sd.models.get_model(image2d)
    assert schema == sd.models.Image2DModel

    # Check the correct data exists.
    result = image2d.data.compute()
    np.testing.assert_equal(result, sample_2d_data[level])

    # Check the metadata.
    metadata = dict(image2d.attrs)
    assert len(metadata) == 1
    assert metadata["transform"] == {"scene0": expected_transformation}