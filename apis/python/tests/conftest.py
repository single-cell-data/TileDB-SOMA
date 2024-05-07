from tempfile import TemporaryDirectory

import anndata
import pytest

import tiledbsoma
import tiledbsoma.io

from ._util import TESTDATA


@pytest.fixture
def conftest_h5ad_path(request):
    # pbmc-small is faster for automated unit-test / CI runs.
    return TESTDATA / "pbmc-small.h5ad"


@pytest.fixture
def conftest_adata(conftest_h5ad_path):
    return anndata.read_h5ad(conftest_h5ad_path)


@pytest.fixture
def conftest_h5ad_file_extended(request):
    # This has more component arrays in it
    return TESTDATA / "pbmc3k_processed.h5ad"


@pytest.fixture
def conftest_adata_extended(conftest_h5ad_file_extended):
    return anndata.read_h5ad(conftest_h5ad_file_extended)


@pytest.fixture
def conftest_pbmc_small(conftest_h5ad_path):
    """Ingest an ``AnnData``, yield a ``TestCase`` with the original and new AnnData objects."""
    with TemporaryDirectory() as exp_path:
        tiledbsoma.io.from_h5ad(exp_path, conftest_h5ad_path, measurement_name="RNA")
        with tiledbsoma.Experiment.open(exp_path) as exp:
            yield exp
