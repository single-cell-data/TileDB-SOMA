import tempfile
from pathlib import Path

import anndata
import pytest

import tiledbsoma
import tiledbsoma.io

HERE = Path(__file__).parent


@pytest.fixture
def h5ad_file(request):
    # pbmc-small is faster for automated unit-test / CI runs.
    # input_path = HERE.parent / "anndata/pbmc3k_processed.h5ad"
    input_path = HERE.parent / "anndata/pbmc-small.h5ad"
    return input_path


@pytest.fixture
def adata(h5ad_file):
    return anndata.read_h5ad(h5ad_file)


def test_platform_config(adata):

    # Set up anndata input path and tiledb-group output path
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name

    # platform_config = {"tiledb": TileDBPlatformConfig(X_capacity=8888))
    # platform_config = {"tiledb": {"X_capacity":8888}}
    # platform_config = {"tiledb": {"capacity":8888}}
    # Hi Paul! Two notes:
    # * In `main-old` we had `X_capacity` separate from capacities for every array. Might be worth
    #   doing here in `main`.
    # * The above three commented-out attempts didn't work for me until I read the code.
    #   This is PEBKAC for sure, however, we might do a bit more aggressive shape-checking
    #   to warn others users, who may be as careless as myself, earlier than later.
    platform_config = {"tiledb": {"create": {"capacity": 8888}}}

    # Ingest
    exp = tiledbsoma.Experiment(output_path)
    tiledbsoma.io.from_anndata(exp, adata, "mRNA", platform_config=platform_config)

    with exp.ms["mRNA"].X["data"]._tiledb_open() as X:
        assert X.schema.capacity == 8888

    # soma.TileDBPlatformConfig(obs_extent=999, X_capacity=8888)
