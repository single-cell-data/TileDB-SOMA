import tempfile
from pathlib import Path

import anndata
import pytest
import tiledb

import tiledbsoma
import tiledbsoma.io
import tiledbsoma.options._tiledb_create_options as tco

HERE = Path(__file__).parent


@pytest.fixture
def h5ad_file(request):
    # pbmc-small is faster for automated unit-test / CI runs.
    # input_path = HERE.parent / "testdata/pbmc3k_processed.h5ad"
    input_path = HERE.parent / "testdata/pbmc-small.h5ad"
    return input_path


@pytest.fixture
def adata(h5ad_file):
    return anndata.read_h5ad(h5ad_file)


def test_platform_config(adata):
    # Set up anndata input path and tiledb-group output path
    with tempfile.TemporaryDirectory() as output_path:
        # Ingest
        tiledbsoma.io.from_anndata(
            output_path,
            adata,
            "RNA",
            platform_config={
                "tiledb": {
                    "create": {
                        "capacity": 8888,
                        "offsets_filters": ["RleFilter", "NoOpFilter"],
                        "dims": {
                            "soma_dim_0": {"tile": 6},
                            # Empty filters for soma_dim_1 overrides the default
                            # dimension zstd level defined below.
                            "soma_dim_1": {"filters": []},
                        },
                        "attrs": {"soma_data": {"filters": ["NoOpFilter"]}},
                        "dataframe_dim_zstd_level": 1,
                        "cell_order": "row-major",
                        "tile_order": "col-major",
                        "sparse_nd_array_dim_zstd_level": 2,
                    }
                }
            },
        )

        with tiledbsoma.Experiment.open(output_path) as exp:
            x_data = exp.ms["RNA"].X["data"]
            x_arr = x_data._handle.reader
            x_sch = x_arr.schema
            assert x_sch.capacity == 8888
            assert x_sch.cell_order == "row-major"
            assert x_sch.tile_order == "col-major"
            assert x_sch.offsets_filters == [tiledb.RleFilter(), tiledb.NoOpFilter()]
            assert x_arr.attr("soma_data").filters == [tiledb.NoOpFilter()]
            assert x_arr.dim("soma_dim_0").tile == 6
            assert x_arr.dim("soma_dim_0").filters == [tiledb.ZstdFilter(level=2)]
            assert x_arr.dim("soma_dim_1").filters == []

            var_df = exp.ms["RNA"].var
            var_arr = var_df._handle.reader
            assert var_arr.dim("soma_joinid").filters == [tiledb.ZstdFilter(level=1)]


def test__from_platform_config__admits_ignored_config_structure():
    try:
        tco.TileDBCreateOptions.from_platform_config(
            dict(
                tiledb=dict(
                    create=tco.TileDBCreateOptions(), future_option=dict(foo="1")
                ),
                not_tiledb=dict(read=dict(buffer_size="128")),
            )
        )
    except Exception as e:
        pytest.fail(f"unexpected exception {e}")


def test__from_platform_config__admits_ignored_options():
    tco.TileDBCreateOptions.from_platform_config(
        {"tiledb": {"create": {"zzz_future_option": "hello"}}}
    )


def test__from_platform_config__admits_plain_dict():
    tdb_create_options = tco.TileDBCreateOptions.from_platform_config(
        {"tiledb": {"create": {"dims": {"soma_dim_0": {"tile": 6}}}}}
    )
    assert tdb_create_options.dim_tile("soma_dim_0") == 6


def test__from_platform_config__admits_create_options_in_dict_shallow():
    tdb_create_options = tco.TileDBCreateOptions.from_platform_config(
        {"tiledb": tco.TileDBCreateOptions(dims={"soma_dim_0": {"tile": 6}})}
    )
    assert tdb_create_options.dim_tile("soma_dim_0") == 6


def test__from_platform_config__admits_create_options_in_dict_at_leaf():
    tdb_create_options = tco.TileDBCreateOptions.from_platform_config(
        {
            "tiledb": {
                "create": tco.TileDBCreateOptions(dims={"soma_dim_0": {"tile": 6}})
            }
        }
    )
    assert tdb_create_options.dim_tile("soma_dim_0") == 6


def test__from_platform_config__admits_create_options_at_root():
    tdb_create_options = tco.TileDBCreateOptions.from_platform_config(
        tco.TileDBCreateOptions(dims={"soma_dim_0": {"tile": 6}})
    )
    assert tdb_create_options.dim_tile("soma_dim_0") == 6


def test_dig_platform_config():
    # Normal.
    assert tco._dig_platform_config(1, int, ("a", "b")) == 1
    assert tco._dig_platform_config({"a": 2}, int, ("a", "b")) == 2
    assert tco._dig_platform_config({"a": {"b": 3}}, int, ("a", "b")) == 3
    assert tco._dig_platform_config(
        {"a": {"b": {"config_data": "hello"}}}, int, ("a", "b")
    ) == {"config_data": "hello"}

    # Missing keys interpolated with empty dict.
    assert tco._dig_platform_config({"x": "y"}, int, ("a", "b")) == {}
    assert tco._dig_platform_config({"a": {"c": "hello"}}, int, ("a", "b")) == {}

    # Unrecognized types (not at tip).
    assert tco._dig_platform_config("hello", int, ("a", "b")) == {}
    assert tco._dig_platform_config({"a": "hello"}, int, ("a", "b")) == {}

    # Unrecognized type (at tip)
    with pytest.raises(TypeError):
        tco._dig_platform_config({"a": {"b": "invalid"}}, int, ("a", "b"))
