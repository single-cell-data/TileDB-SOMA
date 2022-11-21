import tempfile
from pathlib import Path

import anndata
import pytest
import tiledb

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
    with tempfile.TemporaryDirectory() as output_path:
        # Ingest
        exp = tiledbsoma.Experiment(output_path)
        tiledbsoma.io.from_anndata(
            exp,
            adata,
            "mRNA",
            platform_config={
                "tiledb": {
                    "create": {
                        "capacity": 8888,
                        "offsets_filters": ["RleFilter", "NoOpFilter"],
                        "dims": {
                            "soma_dim_0": {"tile": 6},
                            "soma_dim_1": {"filters": []},
                        },
                        "attrs": {"soma_data": {"filters": ["NoOpFilter"]}},
                        "cell_order": "row-major",
                        "tile_order": "col-major",
                    },
                },
            },
        )

        with exp.ms["mRNA"].X["data"]._tiledb_open() as arr:
            sch = arr.schema
            assert sch.capacity == 8888
            assert sch.cell_order == "row-major"
            assert sch.tile_order == "col-major"
            assert sch.offsets_filters == [tiledb.RleFilter(), tiledb.NoOpFilter()]
            assert arr.attr("soma_data").filters == [tiledb.NoOpFilter()]
            assert arr.dim("soma_dim_0").tile == 6
            assert arr.dim("soma_dim_1").filters == []
            print(sch)
