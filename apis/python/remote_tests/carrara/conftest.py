from __future__ import annotations

from collections.abc import Generator
from uuid import uuid4

import anndata as ad
import pytest
import scanpy as sc
import scipy.sparse as sp

import tiledbsoma as soma
import tiledb
import tiledb.client

BASE_URI = "tiledb://TileDB-Inc./Bruce/remote_test"


def pytest_addoption(parser):
    parser.addoption("--carrara", action="store_true", default=False, help="run Carrara tests")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--carrara"):
        # --carrara given in cli: do not skip slow tests
        return
    skip_carrara = pytest.mark.skip(reason="need --carrara option to run")
    for item in items:
        if "carrara" in item.keywords:
            item.add_marker(skip_carrara)


@pytest.fixture(scope="session")
def carrara_context() -> soma.SOMATileDBContext:
    tiledb.client.login(profile_name="qa")
    return soma.SOMATileDBContext(tiledb_ctx=tiledb.Ctx())


@pytest.fixture
def array_path() -> Generator[str, None, None]:
    """Fixture returns an Array path that will be recursively deleted after test finishes."""
    path = f"{BASE_URI}/{uuid4()}"
    yield path

    tiledb.Array.delete_array(path)


@pytest.fixture
def group_path() -> Generator[str, None, None]:
    """Fixture returns a Group path that will be recursively deleted after test finishes."""
    path = f"{BASE_URI}/{uuid4()}"
    yield path

    try:
        with tiledb.Group(path, mode="m") as G:
            G.delete(recursive=True)
    except tiledb.TileDBError:
        pass


@pytest.fixture
def small_pbmc() -> ad.AnnData:
    adata = sc.datasets.pbmc3k_processed()

    # trim it down to something smallish

    adata = adata[0:101, 0:99].copy()
    del adata.obsm["X_draw_graph_fr"]
    del adata.obsm["X_pca"]

    del adata.uns["neighbors"]
    del adata.uns["pca"]
    del adata.uns["rank_genes_groups"]

    del adata.obsp["distances"]

    adata.X = sp.csr_matrix(adata.X)

    return adata
