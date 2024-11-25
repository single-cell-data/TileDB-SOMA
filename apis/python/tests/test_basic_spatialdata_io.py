from urllib.parse import urljoin

import pytest

import tiledbsoma
from tiledbsoma import _factory

spatial_outgest = pytest.importorskip("tiledbsoma.io.spatial.outgest")


def test_outgest_no_spatial(tmp_path, conftest_pbmc_small):
    # Create the SOMA Experiment.
    output_path = urljoin(f"{tmp_path.as_uri()}/", "outgest_no_spatial")
    tiledbsoma.io.from_anndata(output_path, conftest_pbmc_small, measurement_name="RNA")

    # Read full experiment into SpatialData.
    with _factory.open(output_path) as exp:
        sdata = spatial_outgest.to_spatialdata(exp)

    # Check the number of assets (exactly 1 table) is as expected.
    print(sdata)
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
