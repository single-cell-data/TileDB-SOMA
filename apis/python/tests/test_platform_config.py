import json
import tempfile
from pathlib import Path

import pytest

import tiledbsoma
import tiledbsoma.io
import tiledbsoma.options._tiledb_create_write_options as tco

from ._util import assert_adata_equal


def test_platform_config(conftest_pbmc_small):
    # Set up anndata input path and tiledb-group output path
    original = conftest_pbmc_small.copy()
    with tempfile.TemporaryDirectory(prefix="test_platform_config_") as output_path:
        # Ingest
        create_cfg = {
            "capacity": 8888,
            "offsets_filters": [
                "RleFilter",
                {"_type": "GzipFilter", "level": 7},
                "NoOpFilter",
            ],
            "dims": {
                "soma_dim_0": {"tile": 6, "filters": ["RleFilter"]},
                # Empty filters for soma_dim_1 overrides the default
                # dimension zstd level defined below.
                "soma_dim_1": {"filters": []},
            },
            "attrs": {"soma_data": {"filters": ["NoOpFilter"]}},
            "dataframe_dim_zstd_level": 1,
            "cell_order": "row-major",
            "tile_order": "column-major",
            "dense_nd_array_dim_zstd_level": 2,
        }

        tiledbsoma.io.from_anndata(
            output_path,
            conftest_pbmc_small,
            "RNA",
            platform_config={"tiledb": {"create": create_cfg}},
        )
        assert_adata_equal(original, conftest_pbmc_small)

        x_arr_uri = str(Path(output_path) / "ms" / "RNA" / "X" / "data")
        with tiledbsoma.SparseNDArray.open(x_arr_uri) as x_arr:
            cfg = x_arr.schema_config_options()
            assert cfg.capacity == create_cfg["capacity"]
            assert cfg.cell_order == create_cfg["cell_order"]
            assert cfg.tile_order == create_cfg["tile_order"]
            assert json.loads(cfg.offsets_filters) == [
                {"COMPRESSION_LEVEL": -1, "name": "RLE"},
                {"COMPRESSION_LEVEL": 7, "name": "GZIP"},
                {"name": "NOOP"},
            ]

            assert json.loads(cfg.attrs)["soma_data"]["filters"] == [{"name": "NOOP"}]

            soma_dim_0 = json.loads(cfg.dims)["soma_dim_0"]
            assert int(soma_dim_0["tile"]) == 6
            assert soma_dim_0["filters"] == [{"COMPRESSION_LEVEL": -1, "name": "RLE"}]

            # As of 2.17.0 this is the default when empty filter-list, or none at all,
            # is requested. Those who want truly no filtering can request a no-op filter.
            assert json.loads(cfg.dims)["soma_dim_1"]["filters"] == [
                {"COMPRESSION_LEVEL": -1, "name": "ZSTD"}
            ]

        var_arr_uri = str(Path(output_path) / "ms" / "RNA" / "var")
        with tiledbsoma.DataFrame.open(var_arr_uri) as var_arr:
            cfg = var_arr.schema_config_options()
            assert json.loads(cfg.dims)["soma_joinid"]["filters"] == [
                {"COMPRESSION_LEVEL": 1, "name": "ZSTD"}
            ]


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
