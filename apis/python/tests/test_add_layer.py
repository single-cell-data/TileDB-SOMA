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

    soma.X.add_layer_from_matrix_and_dim_values(
        csr, soma.obs.ids(), soma.var.ids(), "data2"
    )

    csr2 = soma.X.data2.csr()
    assert csr2.shape == orig_shape

    tempdir.cleanup()
