import tempfile

import tiledbsoma
import tiledbsoma.io


def test_update_matrix_X(adata_extended):
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name

    tiledbsoma.io.from_anndata(output_path, adata_extended, measurement_name="RNA")

    with tiledbsoma.Experiment.open(output_path) as exp:
        old = exp.ms["RNA"].X["data"].read().tables().concat()

    assert len(old["soma_dim_0"]) == 4848644
    assert len(old["soma_dim_1"]) == 4848644
    assert len(old["soma_data"]) == 4848644

    with tiledbsoma.Experiment.open(output_path, "w") as exp:
        tiledbsoma.io.update_matrix(
            exp.ms["RNA"].X["data"],
            adata_extended.X + 1,
        )

    with tiledbsoma.Experiment.open(output_path) as exp:
        new = exp.ms["RNA"].X["data"].read().tables().concat()

    assert len(new["soma_dim_0"]) == 4848644
    assert len(new["soma_dim_1"]) == 4848644
    assert len(new["soma_data"]) == 4848644

    assert old["soma_dim_0"] == new["soma_dim_0"]
    assert old["soma_dim_1"] == new["soma_dim_1"]
    assert old["soma_data"] != new["soma_data"]


def test_update_matrix_obsm(adata_extended):
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name

    tiledbsoma.io.from_anndata(output_path, adata_extended, measurement_name="RNA")

    with tiledbsoma.Experiment.open(output_path) as exp:
        old = exp.ms["RNA"].obsm["X_pca"].read().tables().concat()

    assert len(old["soma_dim_0"]) == 131900
    assert len(old["soma_dim_1"]) == 131900
    assert len(old["soma_data"]) == 131900

    with tiledbsoma.Experiment.open(output_path, "w") as exp:
        tiledbsoma.io.update_matrix(
            exp.ms["RNA"].obsm["X_pca"],
            adata_extended.obsm["X_pca"] + 1,
        )

    with tiledbsoma.Experiment.open(output_path) as exp:
        new = exp.ms["RNA"].obsm["X_pca"].read().tables().concat()

    assert len(new["soma_dim_0"]) == 131900
    assert len(new["soma_dim_1"]) == 131900
    assert len(new["soma_data"]) == 131900

    assert old["soma_dim_0"] == new["soma_dim_0"]
    assert old["soma_dim_1"] == new["soma_dim_1"]
    assert old["soma_data"] != new["soma_data"]
