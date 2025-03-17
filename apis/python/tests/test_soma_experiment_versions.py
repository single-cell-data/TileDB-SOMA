import os

import pytest

import tiledbsoma

from ._util import ROOT_DATA_DIR


@pytest.mark.parametrize("version", ["1.7.3", "1.12.3", "1.14.5", "1.15.0", "1.15.7"])
@pytest.mark.parametrize(
    "name_and_expected_shape",
    [["pbmc3k_unprocessed", (2700, 13714)], ["pbmc3k_processed", (2638, 1838)]],
)
def test_to_anndata(version, name_and_expected_shape):
    """Checks that experiments written by older versions are still readable,
    in the particular form of doing an outgest."""

    name, expected_shape = name_and_expected_shape
    path = ROOT_DATA_DIR / "soma-experiment-versions" / version / name
    uri = str(path)
    if not os.path.isdir(uri):
        raise RuntimeError(
            f"Missing '{uri}' directory. Try running `make data` "
            "from the TileDB-SOMA project root directory."
        )

    with tiledbsoma.Experiment.open(uri) as exp:
        adata = tiledbsoma.io.to_anndata(
            exp,
            measurement_name="RNA",
            X_layer_name="data",
        )
        assert adata.shape == expected_shape

        expected_nobs = expected_shape[0]
        expected_nvar = expected_shape[1]

        assert adata.obs.shape[0] == expected_nobs

        assert adata.var.shape[0] == expected_nvar

        assert adata.X.shape[0] == expected_nobs
        assert adata.X.shape[1] == expected_nvar

        if name == "pbmc3k_processed":
            for key in ["X_pca"]:
                assert adata.obsm[key].shape == (expected_nobs, 50)

            for key in ["X_tsne", "X_umap", "X_draw_graph_fr"]:
                assert adata.obsm[key].shape == (expected_nobs, 2)

            for key in ["connectivities", "distances"]:
                assert adata.obsp[key].shape == (expected_nobs, expected_nobs)

            assert adata.varm["PCs"].shape == (expected_nvar, 50)
