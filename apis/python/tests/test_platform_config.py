import tempfile

import pytest

import tiledbsoma
import tiledbsoma.io
import tiledbsoma.options._tiledb_create_options as tco
from tiledbsoma._util import verify_obs_var
import tiledb


@pytest.mark.skip(reason="No longer return ArraySchema - see note in test")
def test_platform_config(adata):
    # TODO as we remove usage of TileDB-Py in favor of ArrowSchema, we
    # need a new method to get which filters have applied to the column
    # rather than grabbing it from the ArraySchema. One consideration
    # would be to store TileDB information in JSON format as a field in
    # the ArraySchema metadata very similar to how Pandas stores information
    # within pa.Schema.pandas_metadata. This could hold not only which
    # filters have been applied to the column, but other info that cannot
    # be "directly" stored in the ArrowSchema such as whether the column
    # is a TileDB attribute or dimension, whether this represent a dense
    # or sparse array, etc. This may be as easy as simply copying the
    # platform_config by calling pa.Schema.with_metadata(platform_config).

    # Set up anndata input path and tiledb-group output path
    original = adata.copy()
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
        verify_obs_var(original, adata)

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
            # As of 2.17.0 this is the default when empty filter-list, or none at all,
            # is requested. Those who want truly no filtering can request a no-op filter.
            assert list(x_arr.dim("soma_dim_1").filters) == [
                tiledb.ZstdFilter(level=-1)
            ]
            # TODO as we remove usage of TileDB-Py in favor of ArrowSchema, we
            # need a new method to get which filters have applied to the column
            # rather than grabbing it from the ArraySchema. One consideration
            # would be to store TileDB information in JSON format as a field in
            # the ArraySchema metadata very similar to how Pandas stores information
            # within pa.Schema.pandas_metadata. This could hold not only which
            # filters have been applied to the column, but other info that cannot
            # be "directly" stored in the ArrowSchema such as whether the column
            # is a TileDB attribute or dimension, whether this represent a dense
            # or sparse array, etc. This may be as easy as simply copying the
            # platform_config by calling pa.Schema.with_metadata(platform_config).
            # var_df = exp.ms["RNA"].var
            # var_arr = var_df._handle.reader
            # assert var_arr.dim("soma_joinid").filters == [tiledb.ZstdFilter(level=1)]


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
