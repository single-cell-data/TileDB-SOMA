from tempfile import TemporaryDirectory

import anndata
import pytest

import tiledbsoma
import tiledbsoma.io

from ._util import TESTDATA


@pytest.fixture
def h5ad_file(request):
    # pbmc-small is faster for automated unit-test / CI runs.
    return TESTDATA / "pbmc-small.h5ad"


@pytest.fixture
def adata(h5ad_file):
    return anndata.read_h5ad(h5ad_file)


@pytest.fixture
def h5ad_file_extended(request):
    # This has more component arrays in it
    return TESTDATA / "pbmc3k_processed.h5ad"


@pytest.fixture
def adata_extended(h5ad_file_extended):
    return anndata.read_h5ad(h5ad_file_extended)


@pytest.fixture
def pbmc_small(h5ad_file):
    """Ingest an ``AnnData``, yield a ``TestCase`` with the original and new AnnData objects."""
    with TemporaryDirectory() as exp_path:
        tiledbsoma.io.from_h5ad(exp_path, h5ad_file, measurement_name="RNA")
        with tiledbsoma.Experiment.open(exp_path) as exp:
            yield exp
