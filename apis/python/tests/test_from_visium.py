import os

import pytest

from tiledbsoma.io.spatial import VisiumPaths

from ._util import ROOT_DATA_DIR


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


def test_visium_paths_v2(visium_v2_path):
    """Test ``VisiumPaths`` for Visium v2 in standard structure."""
    visium_paths = VisiumPaths.from_base_folder(visium_v2_path)
    assert os.path.isfile(visium_paths.gene_expression)
    assert os.path.isfile(visium_paths.tissue_positions)
    assert visium_paths.fullres_image is None
    assert os.path.isfile(visium_paths.hires_image)
    assert os.path.isfile(visium_paths.lowres_image)
    assert visium_paths.version == (2, 0, 0)
    assert visium_paths.has_image
    assert visium_paths.major_version == 2
