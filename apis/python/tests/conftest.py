from pathlib import Path
from tempfile import TemporaryDirectory

import anndata
import pytest
from anndata import AnnData

import tiledbsoma
import tiledbsoma.io
from tiledbsoma import Experiment

from ._util import TESTDATA


@pytest.fixture
def conftest_pbmc_small_h5ad_path(request) -> Path:
    """Path to a tiny (80x20) h5ad, useful for unit-test / CI runs."""
    return TESTDATA / "pbmc-small.h5ad"


@pytest.fixture
def conftest_pbmc_small(conftest_pbmc_small_h5ad_path) -> AnnData:
    """Tiny (80x20) AnnData, useful for unit-test / CI runs."""
    return anndata.read_h5ad(conftest_pbmc_small_h5ad_path)


@pytest.fixture
def conftest_pbmc_small_exp(conftest_pbmc_small_h5ad_path) -> Experiment:
    """Ingest an ``AnnData``, yield a ``TestCase`` with the original and new AnnData objects."""
    with TemporaryDirectory() as exp_path:
        tiledbsoma.io.from_h5ad(
            exp_path, conftest_pbmc_small_h5ad_path, measurement_name="RNA"
        )
        with tiledbsoma.Experiment.open(exp_path) as exp:
            yield exp


@pytest.fixture
def conftest_pbmc3k_h5ad_path(request) -> Path:
    """Path to a larger (2638x1838) h5ad, which also includes obsm, obsp, and varm arrays."""
    return TESTDATA / "pbmc3k_processed.h5ad"


@pytest.fixture
def conftest_pbmc3k_adata(conftest_pbmc3k_h5ad_path):
    """Larger (2638x1838) AnnData, which also includes obsm, obsp, and varm arrays."""
    return anndata.read_h5ad(conftest_pbmc3k_h5ad_path)
