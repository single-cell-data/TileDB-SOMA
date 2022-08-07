import tempfile
from pathlib import Path

import anndata
import pytest

import tiledbsc
import tiledbsc.io

HERE = Path(__file__).parent


@pytest.fixture
def h5ad_file(request):
    # Tests in this file rely on specific values form this particular input data file.
    input_path = HERE.parent / "anndata/pbmc-small.h5ad"
    return input_path


@pytest.fixture
def adata(h5ad_file):
    return anndata.read_h5ad(h5ad_file)


def test_or_query(adata):

    # Set up anndata input path and tiledb-group output path
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name

    # Ingest
    soma = tiledbsc.SOMA(output_path)
    tiledbsc.io.from_anndata(soma, adata)

    assert soma.obs.df(attrs=["groups"]).size == 80
    assert soma.obs.query('groups == "g1"', attrs=["groups"]).size == 44
    assert soma.obs.query('groups == "g2"', attrs=["groups"]).size == 36
    assert (
        soma.obs.query('groups == "g1" or groups == "g2"', attrs=["groups"]).size == 80
    )
    assert (
        soma.obs.query('groups == "g1" and groups == "g2"', attrs=["groups"]).size == 0
    )

    tempdir.cleanup()
