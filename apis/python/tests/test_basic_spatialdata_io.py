from urllib.parse import urljoin

import numpy as np
import pyarrow as pa
import pytest

import tiledbsoma as soma
from tiledbsoma import _factory

spatial_outgest = pytest.importorskip("tiledbsoma.io.spatial.outgest")


@pytest.fixture(scope="module")
def experiment_with_single_scene(tmp_path_factory) -> soma.Experiment:
    uri = tmp_path_factory.mktemp("experiment_with_spatial_data").as_uri()
    with soma.Experiment.create(uri) as exp:
        assert exp.uri == uri
        # Create spatial folder.
        exp.add_new_collection("spatial")

        # Create scene 1.
        scene1_uri = urljoin(exp.spatial.uri, "scene1")
        exp.spatial["scene1"] = soma.Scene.create(scene1_uri)
        scene1 = exp.spatial["scene1"]
        assert scene1_uri == scene1.uri
        scene1.coordinate_space = soma.CoordinateSpace.from_axis_names(
            ["x_scene1", "y_scene1"]
        )
        scene1.add_new_collection("obsl")
        scene1.add_new_collection("varl")
        scene1.varl.add_new_collection("RNA")
        scene1.add_new_collection("img")

        # Add point cloud with shape to scene 1.
        points1 = scene1.add_new_point_cloud_dataframe(
            "points1",
            "obsl",
            transform=soma.UniformScaleTransform(
                ("x_scene1", "y_scene1"), ("x", "y"), 2.0
            ),
            schema=pa.schema([("x", pa.float64()), ("y", pa.float64())]),
            domain=[[0, 1], [0, 1]],
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

        # Add point cloud without shape to scene 1.
        points2 = scene1.add_new_point_cloud_dataframe(
            "points2",
            ["varl", "RNA"],
            transform=soma.UniformScaleTransform(
                ("x_scene1", "y_scene1"), ("x", "y"), -1.0
            ),
            schema=pa.schema([("x", pa.float64()), ("y", pa.float64())]),
            domain=[[-1, 0], [-1, 0]],
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

        # Add multiscale image with a single image.
        scene1.add_new_multiscale_image(
            "image1",
            "img",
            type=pa.uint8(),
            level_shape=(3, 64, 64),
            transform=soma.UniformScaleTransform(
                ("x_scene1", "y_scene1"), ("x", "y"), 0.5
            ),
        )
        # TODO: Write data.

        # Add multiscale image with multiple resolutions.
        image2 = scene1.add_new_multiscale_image(
            "image2",
            "img",
            type=pa.uint8(),
            level_key="fullres",
            level_shape=(3, 32, 32),
            transform=soma.UniformScaleTransform(
                ("x_scene1", "y_scene1"), ("x", "y"), 0.5
            ),
        )
        image2.add_new_level("hires", shape=(3, 16, 16))
        image2.add_new_level("lowres", shape=(3, 8, 8))
        scene1.close()

    return soma.Experiment.open(uri, mode="r")


def test_outgest_no_spatial(tmp_path, conftest_pbmc_small):
    # Create the SOMA Experiment.
    output_path = urljoin(f"{tmp_path.as_uri()}/", "outgest_no_spatial")
    soma.io.from_anndata(output_path, conftest_pbmc_small, measurement_name="RNA")

    # Read full experiment into SpatialData.
    with _factory.open(output_path) as exp:
        sdata = spatial_outgest.to_spatial_data(exp)

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


def test_outgest_spatial_only(experiment_with_single_scene):
    # Export to SpatialData.
    sdata = spatial_outgest.to_spatial_data(experiment_with_single_scene)

    # Check the number assets is correct.
    assert len(sdata.tables) == 0
    assert len(sdata.images) == 2
    assert len(sdata.shapes) == 1
    assert len(sdata.points) == 1

    # Check the values of the points.
    assert False  # TODO
