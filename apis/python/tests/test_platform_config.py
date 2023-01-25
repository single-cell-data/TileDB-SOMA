import tempfile
from pathlib import Path

import anndata
import pytest
import tiledb

import tiledbsoma
import tiledbsoma.io
from tiledbsoma.options import TileDBCreateOptions

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
            "RNA",
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
                    }
                }
            },
        )

        with exp.ms["RNA"].X["data"].open_legacy() as data:
            arr = data._tiledb_obj
            sch = arr.schema
            assert sch.capacity == 8888
            assert sch.cell_order == "row-major"
            assert sch.tile_order == "col-major"
            assert sch.offsets_filters == [tiledb.RleFilter(), tiledb.NoOpFilter()]
            assert arr.attr("soma_data").filters == [tiledb.NoOpFilter()]
            assert arr.dim("soma_dim_0").tile == 6
            assert arr.dim("soma_dim_1").filters == []
            print(sch)


def test__from_platform_config__admits_ignored_config_structure():
    try:
        TileDBCreateOptions.from_platform_config(
            dict(
                tiledb=dict(create=TileDBCreateOptions(), future_option=dict(foo="1")),
                not_tiledb=dict(read=dict(buffer_size="128")),
            )
        )
    except Exception as e:
        pytest.fail(f"unexpected exception {e}")


def test__from_platform_config__admits_plain_dict():
    tdb_create_options = TileDBCreateOptions.from_platform_config(
        {"tiledb": {"create": {"dims": {"soma_dim_0": {"tile": 6}}}}}
    )
    assert tdb_create_options.dim_tile("soma_dim_0") == 6


def test__from_platform_config__admits_tiledb_create_options_object():
    tdb_create_options = TileDBCreateOptions.from_platform_config(
        {
            "tiledb": {
                "create": TileDBCreateOptions({"dims": {"soma_dim_0": {"tile": 6}}})
            }
        }
    )
    assert tdb_create_options.dim_tile("soma_dim_0") == 6
