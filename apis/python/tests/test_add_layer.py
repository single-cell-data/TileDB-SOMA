import tempfile
from pathlib import Path

import anndata
import pytest

import tiledbsc
import tiledbsc.io

HERE = Path(__file__).parent


@pytest.fixture
def h5ad_file(request):
    # pbmc-small is faster for automated unit-test / CI runs.
    # input_path = HERE.parent / "anndata/pbmc3k_processed.h5ad"
    input_path = HERE.parent / "anndata/pbmc-small.h5ad"
    return input_path


@pytest.fixture
def adata(h5ad_file):
    return anndata.read_h5ad(h5ad_file)


def test_add_layer(adata):

    # Set up anndata input path and tiledb-group output path
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name

    orig = adata

    # Ingest
    soma = tiledbsc.SOMA(output_path)
    tiledbsc.io.from_anndata(soma, orig)

    csr = soma.X.data.csr()
    orig_shape = csr.shape

    num_rows, num_cols = orig_shape
    for i in range(0, num_rows):
        for j in range(0, num_cols):
            csr[i, j] = i + j

    # Add X layer
    soma.X.add_layer_from_matrix_and_dim_values(
        csr, soma.obs.ids(), soma.var.ids(), "data2"
    )

    csr2 = soma.X.data2.csr()
    assert csr2.shape == orig_shape

    # Add X layer -- more user-friendly syntax
    soma.add_X_layer(csr, "data3")

    csr3 = soma.X.data3.csr()
    assert csr3.shape == orig_shape

    # Add obsm matrix
    soma.obsm.add_matrix_from_matrix_and_dim_values(
        soma.obsm.X_tsne.df(), soma.obs_keys(), "voila"
    )
    assert sorted(soma.obsm.keys()) == ["X_pca", "X_tsne", "voila"]
    assert soma.obsm.voila.shape() == soma.obsm.X_tsne.shape()

    # Add obsp matrix
    soma.obsp.add_matrix_from_matrix_and_dim_values(
        soma.obsp.distances.csr(),
        soma.obs_keys(),
        "voici",
    )
    assert sorted(soma.obsp.keys()) == ["distances", "voici"]
    assert soma.obsp.voici.shape() == soma.obsp.distances.shape()

    tempdir.cleanup()
