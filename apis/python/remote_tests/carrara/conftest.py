from __future__ import annotations

import os
from collections.abc import Generator
from uuid import uuid4

import anndata as ad
import pytest
import scanpy as sc
import scipy.sparse as sp

import tiledbsoma as soma
import tiledb

# Base Carrara URI used for all tests. NB: the teamspace must be owned by the
# test user, allowing access to the S3 path (using tiledb.client API).
#
PROFILE_NAME = os.getenv("CARRARA_TEST_PROFILE") or "qa"
WORKSPACE_NAME = os.getenv("CARRARA_TEST_WORKSPACE") or "TileDB-Inc."
TEAMSPACE_NAME = os.getenv("CARRARA_TEST_TEAMSPACE") or "bruce-uat"
TEST_FOLDER = os.getenv("CARRARA_TEST_FOLDER") or "remote_test"
BASE_URI = f"tiledb://{WORKSPACE_NAME}/{TEAMSPACE_NAME}/{TEST_FOLDER}"


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
def carrara_context() -> soma.SOMAContext:
    import tiledb.client

    tiledb.client.login(profile_name=PROFILE_NAME)
    assert tiledb.client.workspaces.get_workspace(tiledb.client.client.get_workspace_id()).name == WORKSPACE_NAME
    soma.SOMAContext.set_default()
    return soma.SOMAContext.create()


@pytest.fixture
def carrara_array_path() -> Generator[str, None, None]:
    """Fixture returns an Array path that will be recursively deleted after test finishes."""
    import tiledb.client
    path = f"{BASE_URI}/{uuid4()}"
    yield path

    tiledb.client.assets.delete_asset(path, delete_storage=True)


@pytest.fixture
def carrara_group_path() -> Generator[str, None, None]:
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
    adata = adata[0 : adata.n_obs // 2, 0 : adata.n_vars // 2].copy()
    del adata.obsm["X_draw_graph_fr"]
    del adata.obsm["X_pca"]

    del adata.uns["neighbors"]
    del adata.uns["pca"]
    del adata.uns["rank_genes_groups"]

    del adata.obsp["distances"]

    # make the main data matrix sparse
    adata.X = sp.csr_matrix(adata.X)

    # add boolean column to obs for testing
    adata.obs["is_b_cell"] = adata.obs["louvain"] == "B cells"

    return adata
