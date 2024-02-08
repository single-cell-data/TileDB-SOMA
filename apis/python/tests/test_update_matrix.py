import tempfile
from pathlib import Path

import anndata
import pytest

import tiledbsoma
import tiledbsoma.io

HERE = Path(__file__).parent


@pytest.fixture
def h5ad_file(request):
    input_path = HERE.parent / "testdata/pbmc3k_processed.h5ad"
    return input_path


@pytest.fixture
def adata(h5ad_file):
    return anndata.read_h5ad(h5ad_file)


def test_update_matrix_X(adata):
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name

    tiledbsoma.io.from_anndata(output_path, adata, measurement_name="RNA")

    with tiledbsoma.Experiment.open(output_path) as exp:
        old = exp.ms["RNA"].X["data"].read().tables().concat()

    assert len(old["soma_dim_0"]) == 4848644
    assert len(old["soma_dim_1"]) == 4848644
    assert len(old["soma_data"]) == 4848644

    with tiledbsoma.Experiment.open(output_path, "w") as exp:
        tiledbsoma.io.update_matrix(
            exp.ms["RNA"].X["data"],
            adata.X + 1,
        )

    with tiledbsoma.Experiment.open(output_path) as exp:
        new = exp.ms["RNA"].X["data"].read().tables().concat()

    assert len(new["soma_dim_0"]) == 4848644
    assert len(new["soma_dim_1"]) == 4848644
    assert len(new["soma_data"]) == 4848644

    assert old["soma_dim_0"] == new["soma_dim_0"]
    assert old["soma_dim_1"] == new["soma_dim_1"]
    assert old["soma_data"] != new["soma_data"]


def test_update_matrix_obsm(adata):
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name

    tiledbsoma.io.from_anndata(output_path, adata, measurement_name="RNA")

    with tiledbsoma.Experiment.open(output_path) as exp:
        old = exp.ms["RNA"].obsm["X_pca"].read().tables().concat()

    assert len(old["soma_dim_0"]) == 131900
    assert len(old["soma_dim_1"]) == 131900
    assert len(old["soma_data"]) == 131900

    with tiledbsoma.Experiment.open(output_path, "w") as exp:
        tiledbsoma.io.update_matrix(
            exp.ms["RNA"].obsm["X_pca"],
            adata.obsm["X_pca"] + 1,
        )

    with tiledbsoma.Experiment.open(output_path) as exp:
        new = exp.ms["RNA"].obsm["X_pca"].read().tables().concat()

    assert len(new["soma_dim_0"]) == 131900
    assert len(new["soma_dim_1"]) == 131900
    assert len(new["soma_data"]) == 131900

    assert old["soma_dim_0"] == new["soma_dim_0"]
    assert old["soma_dim_1"] == new["soma_dim_1"]
    assert old["soma_data"] != new["soma_data"]
