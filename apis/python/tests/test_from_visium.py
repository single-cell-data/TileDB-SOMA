import os

import numpy as np
import pandas as pd
import pytest

import tiledbsoma as soma

from ._util import ROOT_DATA_DIR

spatial_io = pytest.importorskip("tiledbsoma.io.spatial")


@pytest.fixture(scope="module")
def visium_v1_path():
    """Fixture that checks the example Visium v1 dataset exists."""
    visium_path = ROOT_DATA_DIR / "example-visium-v1"
    if not os.path.isdir(visium_path):
        raise RuntimeError(
            "Missing 'data/example-visium-v1' directory. Try running `make data` "
            "from the TileDB-SOMA project root directory."
        )
    for filename in [
        "filtered_feature_bc_matrix.h5",
        "raw_feature_bc_matrix.h5",
        "spatial/tissue_positions_list.csv",
        "filtered_visium1_loc.csv",
        "spatial/scalefactors_json.json",
        "spatial/tissue_hires_image.png",
        "spatial/tissue_lowres_image.png",
    ]:
        if not os.path.isfile(visium_path / filename):
            raise RuntimeError(
                f"Missing file 'data/example-visium-v1/{filename}'. Try removing "
                f"the directory 'data/example-visium-v1' and re-running `make data'"
                f"from the project root directory."
            )

    return visium_path


@pytest.fixture(scope="module")
def visium_v2_path():
    """Fixture that checks the example Visium v2 dataset exists."""
    visium_path = ROOT_DATA_DIR / "example-visium-v2"
    if not os.path.isdir(visium_path):
        raise RuntimeError(
            "Missing 'data/example-visium-v2' directory. Try running `make data` "
            "from the TileDB-SOMA project root directory."
        )
    for filename in [
        "filtered_feature_bc_matrix.h5",
        "raw_feature_bc_matrix.h5",
        "filtered_visium2_loc.csv",
        "spatial/tissue_positions.csv",
        "spatial/scalefactors_json.json",
        "spatial/tissue_hires_image.png",
        "spatial/tissue_lowres_image.png",
    ]:
        if not os.path.isfile(visium_path / filename):
            raise RuntimeError(
                f"Missing file 'data/example-visium-v2/{filename}'. Try removing "
                f"the directory 'data/example-visium-v2' and re-running `make data'"
                f"from the project root directory."
            )

    return visium_path


def test_visium_paths_v1(visium_v1_path):
    """Test ``VisiumPaths`` for Visium v1 in standard structure."""
    visium_paths = spatial_io.VisiumPaths.from_base_folder(visium_v1_path)
    assert os.path.isfile(visium_paths.gene_expression)
    assert os.path.isfile(visium_paths.tissue_positions)
    assert visium_paths.fullres_image is None
    assert os.path.isfile(visium_paths.hires_image)
    assert os.path.isfile(visium_paths.lowres_image)
    assert visium_paths.version == (1, 1, 0)
    assert visium_paths.has_image
    assert visium_paths.major_version == 1


def test_visium_paths_v2(visium_v2_path):
    """Test ``VisiumPaths`` for Visium v2 in standard structure."""
    visium_paths = spatial_io.VisiumPaths.from_base_folder(visium_v2_path)
    assert os.path.isfile(visium_paths.gene_expression)
    assert os.path.isfile(visium_paths.tissue_positions)
    assert visium_paths.fullres_image is None
    assert os.path.isfile(visium_paths.hires_image)
    assert os.path.isfile(visium_paths.lowres_image)
    assert visium_paths.version == (2, 0, 0)
    assert visium_paths.has_image
    assert visium_paths.major_version == 2


@pytest.mark.slow
def test_from_visium_for_visium_v1(tmp_path, visium_v1_path):
    """Test `from_visium` runs without error."""
    PIL = pytest.importorskip("PIL")
    uri = f"{tmp_path.as_uri()}/from_visium_for_visium_v2"
    exp_uri = spatial_io.from_visium(
        uri,
        visium_v1_path,
        "RNA",
        "human_lymph_node",
        write_obs_spatial_presence=True,
        write_var_spatial_presence=True,
    )
    with soma.Experiment.open(exp_uri) as exp:

        # Check for the existance of obs, RNA/X, and RNA/var
        assert isinstance(exp.obs, soma.DataFrame)
        assert isinstance(exp.ms["RNA"].X["data"], soma.SparseNDArray)
        assert isinstance(exp.ms["RNA"].var, soma.DataFrame)

        # Check for the existance of the presence matrices.
        assert isinstance(exp.obs_spatial_presence, soma.DataFrame)
        assert isinstance(exp.ms["RNA"].var_spatial_presence, soma.DataFrame)

        # Define some expected data.
        pixel_coord_space = soma.CoordinateSpace(
            (soma.Axis("x", "pixels"), soma.Axis("y", "pixels"))
        )

        # Check the scene exists and has the desired metadata.
        assert isinstance(exp.spatial["human_lymph_node"], soma.Scene)
        scene = exp.spatial["human_lymph_node"]
        assert scene.coordinate_space == pixel_coord_space

        # Check the scene subcollections:
        # - `obsl` only has a point cloud dataframe named `loc`,
        # - `varl` is empty,
        # - `img` only has a multiscale image named `tissue`.
        assert len(scene.obsl.items()) == 1
        assert isinstance(scene.obsl["loc"], soma.PointCloudDataFrame)
        assert len(scene.varl.items()) == 0
        assert len(scene.img.items()) == 1
        assert isinstance(scene.img["tissue"], soma.MultiscaleImage)

        # Check transform to `loc`.
        loc_transform = scene.get_transform_to_point_cloud_dataframe("loc", "obsl")
        assert isinstance(loc_transform, soma.IdentityTransform)
        assert loc_transform.input_axes == ("x", "y")
        assert loc_transform.output_axes == ("x", "y")

        # Check transform to `img`.
        img_transform = scene.get_transform_to_multiscale_image("tissue", "img")
        assert isinstance(img_transform, soma.ScaleTransform)
        assert img_transform.input_axes == ("x", "y")
        assert img_transform.output_axes == ("x", "y")
        np.testing.assert_allclose(
            img_transform.scale_factors,
            np.array([0.170105, 0.170111], dtype=np.float64),
            atol=0.000001,
        )

        # Check spot locations.
        loc = scene.obsl["loc"]
        assert loc.coordinate_space == pixel_coord_space
        actual_loc_data = loc.read().concat().to_pandas()
        expected_loc_data = pd.read_csv(visium_v1_path / "filtered_visium1_loc.csv")
        pd.testing.assert_frame_equal(actual_loc_data, expected_loc_data)

        # Check tissue image.
        image = scene.img["tissue"]
        assert image.coordinate_space == pixel_coord_space
        hires_data = np.moveaxis(image["hires"].read().to_numpy(), 0, -1)
        with PIL.Image.open(
            visium_v1_path / "spatial" / "tissue_hires_image.png"
        ) as input_hires:
            expected = np.array(input_hires)
            np.testing.assert_equal(expected, hires_data)
        lowres_data = np.moveaxis(image["lowres"].read().to_numpy(), 0, -1)
        with PIL.Image.open(
            visium_v1_path / "spatial" / "tissue_lowres_image.png"
        ) as input_lowres:
            expected = np.array(input_lowres)
            np.testing.assert_equal(expected, lowres_data)


@pytest.mark.slow
def test_from_visium_for_visium_v2(tmp_path, visium_v2_path):
    """Test `from_visium` runs without error."""
    PIL = pytest.importorskip("PIL")
    uri = f"{tmp_path.as_uri()}/from_visium_for_visium_v2"
    exp_uri = spatial_io.from_visium(
        uri,
        visium_v2_path,
        "RNA",
        "fresh_frozen_mouse_brain",
        write_obs_spatial_presence=True,
        write_var_spatial_presence=True,
    )
    with soma.Experiment.open(exp_uri) as exp:

        # Check for the existance of obs, RNA/X, and RNA/var
        assert isinstance(exp.obs, soma.DataFrame)
        assert isinstance(exp.ms["RNA"].X["data"], soma.SparseNDArray)
        assert isinstance(exp.ms["RNA"].var, soma.DataFrame)

        # Check for the existance of the presence matrices.
        assert isinstance(exp.obs_spatial_presence, soma.DataFrame)
        assert isinstance(exp.ms["RNA"].var_spatial_presence, soma.DataFrame)

        # Define some expected data.
        pixel_coord_space = soma.CoordinateSpace(
            (soma.Axis("x", "pixels"), soma.Axis("y", "pixels"))
        )

        # Check the scene exists and has the desired metadata.
        assert isinstance(exp.spatial["fresh_frozen_mouse_brain"], soma.Scene)
        scene = exp.spatial["fresh_frozen_mouse_brain"]
        assert scene.coordinate_space == pixel_coord_space

        # Check the scene subcollections:
        # - `obsl` only has a point cloud dataframe named `loc`,
        # - `varl` is empty,
        # - `img` only has a multiscale image named `tissue`.
        assert len(scene.obsl.items()) == 1
        assert isinstance(scene.obsl["loc"], soma.PointCloudDataFrame)
        assert len(scene.varl.items()) == 0
        assert len(scene.img.items()) == 1
        assert isinstance(scene.img["tissue"], soma.MultiscaleImage)

        # Check transform to `loc`.
        loc_transform = scene.get_transform_to_point_cloud_dataframe("loc", "obsl")
        assert isinstance(loc_transform, soma.IdentityTransform)
        assert loc_transform.input_axes == ("x", "y")
        assert loc_transform.output_axes == ("x", "y")

        # Check transform to `img`.
        img_transform = scene.get_transform_to_multiscale_image("tissue", "img")
        assert isinstance(img_transform, soma.ScaleTransform)
        assert img_transform.input_axes == ("x", "y")
        assert img_transform.output_axes == ("x", "y")
        np.testing.assert_allclose(
            img_transform.scale_factors,
            np.array([0.150015, 0.150015], dtype=np.float64),
            atol=0.000001,
        )

        # Check point cloud dataframe data.
        loc = scene.obsl["loc"]
        assert loc.coordinate_space == pixel_coord_space
        actual_loc_data = loc.read().concat().to_pandas()
        expected_loc_data = pd.read_csv(visium_v2_path / "filtered_visium2_loc.csv")
        pd.testing.assert_frame_equal(actual_loc_data, expected_loc_data)

        # Check image.
        image = scene.img["tissue"]
        assert image.coordinate_space == pixel_coord_space
        hires_data = np.moveaxis(image["hires"].read().to_numpy(), 0, -1)
        with PIL.Image.open(
            visium_v2_path / "spatial" / "tissue_hires_image.png"
        ) as input_hires:
            expected = np.array(input_hires)
            np.testing.assert_equal(expected, hires_data)
        lowres_data = np.moveaxis(image["lowres"].read().to_numpy(), 0, -1)
        with PIL.Image.open(
            visium_v2_path / "spatial" / "tissue_lowres_image.png"
        ) as input_lowres:
            expected = np.array(input_lowres)
            np.testing.assert_equal(expected, lowres_data)
