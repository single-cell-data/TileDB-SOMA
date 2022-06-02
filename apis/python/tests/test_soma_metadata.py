import anndata
import tiledb
import tiledbsc
import tiledbsc.io

import numpy as np

import pytest
import tempfile
import os
from pathlib import Path

HERE = Path(__file__).parent


@pytest.fixture
def h5ad_file(request):
    input_path = HERE.parent / "anndata/pbmc-small.h5ad"
    return input_path


def test_soma_metadata(h5ad_file):
    """
    Verify basic metadata access at the tiledbsc-py level.
    """

    # Set up anndata input path and tiledb-group output path
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name

    # Ingest
    soma = tiledbsc.SOMA(output_path, verbose=False)
    tiledbsc.io.from_h5ad(soma, h5ad_file)
    assert soma.exists()

    # Group-level metadata
    assert soma.has_metadata("nonesuch") == False

    soma.set_metadata("foo", "bar")
    assert soma.has_metadata("nonesuch") == False
    assert soma.has_metadata("foo") == True
    assert soma.get_metadata("foo") == "bar"

    # Array-level metadata
    assert soma.obs.has_metadata("nonesuch") == False

    soma.obs.set_metadata("foo", "bar")
    assert soma.obs.has_metadata("nonesuch") == False
    assert soma.obs.has_metadata("foo") == True
    assert "foo" in soma.obs.metadata_keys()
    assert soma.obs.get_metadata("foo") == "bar"

    tempdir.cleanup()
