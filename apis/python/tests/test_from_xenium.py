import os

import pytest

import tiledbsoma as soma

from ._util import ROOT_DATA_DIR

spatial_io = pytest.importorskip("tiledbsoma.io.spatial")


@pytest.fixture(scope="module")
def xenium_v1_path():
    """Fixture that checks the example Xenium v1 dataset exists."""
    xenium_path = ROOT_DATA_DIR / "example-xenium-v1"
    if not os.path.isdir(xenium_path):
        raise RuntimeError(
            "Missing 'data/example-xenium-v1' directory. Try running `make data` "
            "from the TileDB-SOMA project root directory."
        )
    for filename in [
        "cell_feature_matrix.h5",
        "cells.parquet",
        "cell_boundaries.parquet",
        "nucleus_boundaries.parquet",
        "experiment.xenium",
    ]:
        if not os.path.isfile(xenium_path / filename):
            raise RuntimeError(
                f"Missing file 'data/example-xenium-v1/{filename}'. Try removing "
                f"the directory 'data/example-xenium-v1' and re-running `make data'"
                f"from the project root directory."
            )

    return xenium_path


@pytest.mark.slow
def test_from_xenium_for_xenium_v1(tmp_path, xenium_v1_path):
    """Test `from_xenium` runs without error."""
    uri = f"{tmp_path.as_uri()}/from_visium_for_xenium_v1"
    exp_uri = spatial_io.from_xenium(
        uri,
        xenium_v1_path,
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
