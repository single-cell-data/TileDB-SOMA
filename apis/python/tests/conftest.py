from __future__ import annotations

import multiprocessing
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

import anndata
import pytest
from anndata import AnnData

import tiledbsoma
import tiledbsoma.io
from tiledbsoma import Experiment

from ._util import TESTDATA


@pytest.fixture(scope="session")
def soma_tiledb_config() -> dict[str, Any] | None:
    # Configuration primarily focuses on memory usage, as the CI environment
    # is often memory constrained. The smallest runners are currently 7GiB
    # of RAM, whereas the TileDB core has a default memory budget exceeding
    # 10GiB.

    tiledb_config: dict | None = None

    is_CI = os.getenv("CI", "false") == "true"
    if is_CI:
        tiledb_config = {
            "sm.mem.total_budget": 1 * 1024**3,
            "sm.memory_budget": 512 * 1024**2,
            "sm.memory_budget_var": 512 * 1024**2,
            "soma.init_buffer_bytes": 128 * 1024**2,
        }
    return tiledb_config


@pytest.fixture(scope="module")
def soma_tiledb_context(soma_tiledb_config: dict[str, Any] | None) -> tiledbsoma.SOMATileDBContext:
    """Fixture which builds a SOMATileDBContext based on a reasonable default configuration."""
    return tiledbsoma.SOMATileDBContext(tiledb_config=soma_tiledb_config)


@pytest.fixture
def conftest_pbmc_small_h5ad_path(request) -> Path:
    """Path to a tiny (80x20) h5ad, useful for unit-test / CI runs."""
    return TESTDATA / "pbmc-small.h5ad"


@pytest.fixture
def conftest_pbmc_small(conftest_pbmc_small_h5ad_path) -> AnnData:
    """Tiny (80x20) AnnData, useful for unit-test / CI runs."""
    return anndata.read_h5ad(conftest_pbmc_small_h5ad_path)


@pytest.fixture
def conftest_pbmc_small_exp_path(conftest_pbmc_small_h5ad_path) -> Path:
    with TemporaryDirectory("conftest_pbmc_small_exp_") as exp_path:
        tiledbsoma.io.from_h5ad(
            exp_path,
            conftest_pbmc_small_h5ad_path,
            measurement_name="RNA",
        )
        yield exp_path


@pytest.fixture
def conftest_pbmc_small_exp(conftest_pbmc_small_h5ad_path: Path) -> Experiment:
    """Ingest an ``AnnData``, yield a ``TestCase`` with the original and new AnnData objects."""
    with TemporaryDirectory("conftest_pbmc_small_exp_") as exp_path:
        tiledbsoma.io.from_h5ad(
            exp_path,
            conftest_pbmc_small_h5ad_path,
            measurement_name="RNA",
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


multiprocessing.set_start_method("spawn", force=True)
