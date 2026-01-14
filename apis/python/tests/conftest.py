from __future__ import annotations

import multiprocessing
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

import anndata as ad
import pytest
from anndata import AnnData

import tiledbsoma
import tiledbsoma.io
from tiledbsoma import Experiment

from ._util import TESTDATA


@pytest.fixture(scope="session")
def soma_tiledb_config() -> dict[str, Any] | None:
    # Configuration primarily focuses on memory usage, as the CI environment
    # is memory constrained. The smallest runners have 7GiB of RAM, whereas
    # whereas the TileDB core has a default memory budget of 10GiB.

    # See https://docs.github.com/en/actions/how-tos/writing-workflows/choosing-what-your-workflow-does/store-information-in-variables
    # for default variables in GHA.

    tiledb_config: dict | None = None

    is_CI = os.getenv("CI", "false") == "true"
    if is_CI:
        # default concurrency is cpu_count. Halve to reduce per-worker memory use
        n_cpus = max(1, (os.cpu_count() or 1) // 4)
        tiledb_config = {
            "sm.mem.total_budget": 512 * 1024**2,
            "sm.memory_budget": 256 * 1024**2,
            "sm.memory_budget_var": 256 * 1024**2,
            "soma.init_buffer_bytes": 64 * 1024**2,
            "sm.compute_concurrency_level": n_cpus,
            "sm.io_concurrency_level": n_cpus,
        }
    return tiledb_config


@pytest.fixture(scope="module")
def soma_tiledb_context(soma_tiledb_config: dict[str, Any] | None) -> tiledbsoma.SOMAContext:
    """Fixture which builds a SOMAContext based on a reasonable (i.e., small) default configuration."""
    return tiledbsoma.SOMAContext.create(config=soma_tiledb_config)


@pytest.fixture
def conftest_pbmc_small_h5ad_path(request) -> Path:
    """Path to a tiny (80x20) h5ad, useful for unit-test / CI runs."""
    if not TESTDATA.exists():
        raise RuntimeError(f"Missing directory '{TESTDATA}'. Try re-running `make data` from the project root.")
    return TESTDATA / "pbmc-small.h5ad"


@pytest.fixture
def conftest_pbmc_small(conftest_pbmc_small_h5ad_path) -> AnnData:
    """Tiny (80x20) AnnData, useful for unit-test / CI runs."""
    return ad.read_h5ad(conftest_pbmc_small_h5ad_path)


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
    if not TESTDATA.exists():
        raise RuntimeError(f"Missing directory '{TESTDATA}'. Try re-running `make data` from the project root.")
    return TESTDATA / "pbmc3k_processed.h5ad"


@pytest.fixture
def conftest_pbmc3k_adata(conftest_pbmc3k_h5ad_path):
    """Larger (2638x1838) AnnData, which also includes obsm, obsp, and varm arrays."""
    return ad.read_h5ad(conftest_pbmc3k_h5ad_path)


multiprocessing.set_start_method("spawn", force=True)
