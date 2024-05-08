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
def pbmc0_h5ad_path(request) -> Path:
    """Path to a tiny (80x20) h5ad, useful for unit-test / CI runs."""
    return TESTDATA / "pbmc-small.h5ad"


@pytest.fixture
def pbmc0_adata(pbmc0_h5ad_path) -> AnnData:
    """Tiny (80x20) AnnData, useful for unit-test / CI runs."""
    return anndata.read_h5ad(pbmc0_h5ad_path)


@pytest.fixture
def pbmc0_exp(pbmc0_h5ad_path) -> Experiment:
    """Ingest an ``AnnData``, yield a ``TestCase`` with the original and new AnnData objects."""
    with TemporaryDirectory() as exp_path:
        tiledbsoma.io.from_h5ad(exp_path, pbmc0_h5ad_path, measurement_name="RNA")
        with tiledbsoma.Experiment.open(exp_path) as exp:
            yield exp


@pytest.fixture
def pbmc_3k_h5ad_path(request) -> Path:
    """Path to a larger (2638x1838) h5ad, which also includes obsm, obsp, and varm arrays."""
    return TESTDATA / "pbmc3k_processed.h5ad"


@pytest.fixture
def pbmc_3k_adata(pbmc_3k_h5ad_path):
    """Larger (2638x1838) AnnData, which also includes obsm, obsp, and varm arrays."""
    return anndata.read_h5ad(pbmc_3k_h5ad_path)
