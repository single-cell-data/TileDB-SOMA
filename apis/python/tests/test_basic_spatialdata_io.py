from urllib.parse import urljoin

import numpy as np
import pandas as pd
import pyarrow as pa
import pytest

import tiledbsoma as soma
from tiledbsoma import _factory

spatial_outgest = pytest.importorskip("tiledbsoma.io.spatial.outgest")
sd = pytest.importorskip("spatialdata")
gpd = pytest.importorskip("geopandas")
shapely = pytest.importorskip("shapely")


@pytest.fixture(scope="module")
def sample_2d_data():
    return [
        np.random.randint(0, 255, size=(3, 64, 64), dtype=np.uint8),
        np.random.randint(0, 255, size=(3, 32, 32), dtype=np.uint8),
        np.random.randint(0, 255, size=(3, 16, 16), dtype=np.uint8),
        np.random.randint(0, 255, size=(3, 8, 8), dtype=np.uint8),
    ]


@pytest.fixture(scope="module")
def experiment_with_single_scene(tmp_path_factory, sample_2d_data) -> soma.Experiment:
    uri = tmp_path_factory.mktemp("experiment_with_spatialdata").as_uri()
    with soma.Experiment.create(uri) as exp:
        assert exp.uri == uri
        # Create spatial folder.
        with exp.add_new_collection("spatial") as spatial:

            # Create scene 1.
            with spatial.add_new_collection(
                "scene1", soma.Scene, coordinate_space=("x_scene1", "y_scene1")
            ) as scene1:
                scene1.add_new_collection("obsl")
                scene1.add_new_collection("varl")
                scene1.varl.add_new_collection("RNA")
                scene1.add_new_collection("img")

                # Add point cloud with shape to scene 1 obsl.
                points1 = scene1.add_new_point_cloud_dataframe(
                    "points1",
                    "obsl",
                    transform=soma.UniformScaleTransform(
                        ("x_scene1", "y_scene1"), ("x", "y"), 2.0
                    ),
                    schema=pa.schema([("x", pa.float64()), ("y", pa.float64())]),
                    domain=[[0, 1], [0, 1], [0, 3]],
                )
                points1.write(
                    pa.Table.from_pydict(
                        {
                            "soma_joinid": np.arange(4),
                            "x": np.array([0, 0, 0.5, 0.5]),
                            "y": np.array([0, 0.5, 0, 0.5]),
                        }
                    )
                )
                points1.metadata["soma_geometry"] = 1.0
                points1.metadata["soma_geometry_type"] = "radius"

                # Add point cloud wihtout shape to scene 1 obsl
                points3 = scene1.add_new_point_cloud_dataframe(
                    "points3",
                    "obsl",
                    transform=soma.UniformScaleTransform(
                        ("x_scene1", "y_scene1"), ("x", "y"), 4.0
                    ),
                    schema=pa.schema([("x", pa.float64()), ("y", pa.float64())]),
                    domain=[[-1, 0], [-1, 0], [0, 3]],
                )
                points3.write(
                    pa.Table.from_pydict(
                        {
                            "soma_joinid": np.arange(4),
                            "x": np.array([0, 0, -0.5, -0.5]),
                            "y": np.array([0, -0.5, 0, -0.5]),
                        }
                    )
                )

                # Add point cloud without shape to scene 1 varl.
                points2 = scene1.add_new_point_cloud_dataframe(
                    "points2",
                    ["varl", "RNA"],
                    transform=soma.UniformScaleTransform(
                        ("x_scene1", "y_scene1"), ("x", "y"), -1.0
                    ),
                    schema=pa.schema([("x", pa.float64()), ("y", pa.float64())]),
                    domain=[[-1, 0], [-1, 0], [0, 3]],
                )
                points2.write(
                    pa.Table.from_pydict(
                        {
                            "soma_joinid": np.arange(4),
                            "x": np.array([0, 0, -0.5, -0.5]),
                            "y": np.array([0, -0.5, 0, -0.5]),
                        }
                    )
                )

                # Add point cloud with shape to scene 1 varl.
                points4 = scene1.add_new_point_cloud_dataframe(
                    "points4",
                    ["varl", "RNA"],
                    transform=soma.UniformScaleTransform(
                        ("x_scene1", "y_scene1"), ("x", "y"), 0.25
                    ),
                    schema=pa.schema([("x", pa.float64()), ("y", pa.float64())]),
                    domain=[[0, 1], [0, 1], [0, 3]],
                )
                points4.write(
                    pa.Table.from_pydict(
                        {
                            "soma_joinid": np.arange(4),
                            "x": np.array([0, 0, 0.5, 0.5]),
                            "y": np.array([0, 0.5, 0, 0.5]),
                        }
                    )
                )
                points4.metadata["soma_geometry"] = 2.0
                points4.metadata["soma_geometry_type"] = "radius"

                # Add multiscale image with a single image.
                with scene1.add_new_multiscale_image(
                    "image1",
                    "img",
                    type=pa.uint8(),
                    level_shape=(3, 64, 64),
                    transform=soma.UniformScaleTransform(
                        ("x_scene1", "y_scene1"), ("x", "y"), 0.5
                    ),
                ) as image1:
                    coords = (slice(None), slice(None), slice(None))
                    l0 = image1["level0"]
                    l0.write(coords, pa.Tensor.from_numpy(sample_2d_data[0]))

                # Add multiscale image with multiple resolutions.
                with scene1.add_new_multiscale_image(
                    "image2",
                    "img",
                    type=pa.uint8(),
                    level_key="fullres",
                    level_shape=(3, 32, 32),
                    transform=soma.UniformScaleTransform(
                        ("x_scene1", "y_scene1"), ("x", "y"), 0.5
                    ),
                ) as image2:
                    coords = (slice(None), slice(None), slice(None))
                    fullres = image2["fullres"]
                    fullres.write(coords, pa.Tensor.from_numpy(sample_2d_data[1]))
                    hires = image2.add_new_level("hires", shape=(3, 16, 16))
                    hires.write(coords, pa.Tensor.from_numpy(sample_2d_data[2]))
                    lowres = image2.add_new_level("lowres", shape=(3, 8, 8))
                    lowres.write(coords, pa.Tensor.from_numpy(sample_2d_data[3]))

    return soma.Experiment.open(uri, mode="r")


def test_outgest_no_spatial(tmp_path, conftest_pbmc_small):
    # Create the SOMA Experiment.
    output_path = urljoin(f"{tmp_path.as_uri()}/", "outgest_no_spatial")
    soma.io.from_anndata(output_path, conftest_pbmc_small, measurement_name="RNA")

    # Read full experiment into SpatialData.
    with _factory.open(output_path) as exp:
        sdata = spatial_outgest.to_spatialdata(exp)

    # Check the number of assets (exactly 1 table) is as expected.
    assert len(sdata.tables) == 2
    assert len(sdata.points) == 0
    assert len(sdata.shapes) == 0
    assert len(sdata.images) == 0

    # Check the values of the anndata table.
    rna = sdata.tables["RNA"]
    assert rna.obs.shape == conftest_pbmc_small.obs.shape
    assert rna.var.shape == conftest_pbmc_small.var.shape
    assert rna.X.shape == conftest_pbmc_small.X.shape

    for key in conftest_pbmc_small.obsm.keys():
        assert rna.obsm[key].shape == conftest_pbmc_small.obsm[key].shape
    for key in conftest_pbmc_small.varm.keys():
        assert rna.varm[key].shape == conftest_pbmc_small.varm[key].shape
    for key in conftest_pbmc_small.obsp.keys():
        assert rna.obsp[key].shape == conftest_pbmc_small.obsp[key].shape
    for key in conftest_pbmc_small.varp.keys():
        assert rna.varp[key].shape == conftest_pbmc_small.varp[key].shape

    # Check the values of the anndata table.
    raw = sdata.tables["raw"]
    assert raw.var.shape == conftest_pbmc_small.raw.var.shape
    assert raw.X.shape == conftest_pbmc_small.raw.shape


def test_outgest_spatial_only(experiment_with_single_scene, sample_2d_data):
    # Export to SpatialData.
    sdata = spatial_outgest.to_spatialdata(experiment_with_single_scene)

    # Check the number of assets is correct.
    assert len(sdata.tables) == 0
    assert len(sdata.images) == 2
    assert len(sdata.shapes) == 2
    assert len(sdata.points) == 2

    # Check image1.
    image1 = sdata.images["scene1_image1"]
    assert isinstance(image1, sd.models.models.DataArray)  # Verify single scale image.
    assert image1.attrs["transform"] == {
        "scene1": sd.transformations.Scale([2, 2], ("x", "y"))
    }
    image_data = image1.data.compute()
    np.testing.assert_equal(image_data, sample_2d_data[0])

    # Check image2.
    image2 = sdata.images["scene1_image2"]
    assert isinstance(image2, sd.models.models.DataTree)  # Verify mulitscale image.
    for index in range(3):
        image_level = image2[f"scale{index}"]["image"]
        image_data = image_level.data.compute()
        np.testing.assert_equal(image_data, sample_2d_data[index + 1])
        scale = 2 ** (index + 1)
        assert image_level.attrs["transform"] == {
            "scene1": sd.transformations.Scale([scale, scale], ("x", "y"))
        }

    # Check points1.
    points1 = sdata.shapes["scene1_points1"]
    points1_expected = gpd.GeoDataFrame.from_dict(
        {
            "obs_id": np.arange(4),
            "radius": np.ones((4,), dtype=np.float64),
            "geometry": shapely.points(
                [[0, 0], [0, 0.5], [0.5, 0], [0.5, 0.5]]
            ).tolist(),
        }
    )
    assert all(points1 == points1_expected)
    assert points1.attrs["transform"] == {
        "scene1": sd.transformations.Scale([0.5, 0.5], ("x", "y"))
    }

    # Check points2.
    points2 = sdata.points["scene1_RNA_points2"]
    points2_expected = pd.DataFrame.from_dict(
        {
            "x": np.array([0, 0, -0.5, -0.5]),
            "y": np.array([0, -0.5, 0, -0.5]),
            "var_id": np.arange(4),
        }
    )
    points2_data = points2.compute()
    print(points2_data)
    assert all(points2_data == points2_expected)
    assert points2.attrs["transform"] == {
        "scene1": sd.transformations.Scale([-1.0, -1.0], ("x", "y"))
    }

    # Check points3.
    points3 = sdata.points["scene1_points3"]
    points3_expected = pd.DataFrame.from_dict(
        {
            "x": np.array([0, 0, -0.5, -0.5]),
            "y": np.array([0, -0.5, 0, -0.5]),
            "obs_id": np.arange(4),
        }
    )
    points3_data = points3.compute()
    print(points3_data)
    assert all(points3_data == points3_expected)
    assert points3.attrs["transform"] == {
        "scene1": sd.transformations.Scale([0.25, 0.25], ("x", "y"))
    }

    # Check points4.
    points4 = sdata.shapes["scene1_RNA_points4"]
    points4_expected = gpd.GeoDataFrame.from_dict(
        {
            "var_id": np.arange(4),
            "radius": np.ones((4,), dtype=np.float64),
            "geometry": shapely.points(
                [[0, 0], [0, 0.5], [0.5, 0], [0.5, 0.5]]
            ).tolist(),
        }
    )
    assert all(points4 == points4_expected)
    assert points4.attrs["transform"] == {
        "scene1": sd.transformations.Scale([4.0, 4.0], ("x", "y"))
    }
