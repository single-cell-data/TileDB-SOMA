import numpy as np
import pandas as pd
import pyarrow as pa
import pytest

import tiledbsoma as soma

"""
Test read/write of unicode, ascii and binary
"""


@pytest.fixture
def sample_arrow_table():
    df = pd.DataFrame(
        data={
            "soma_joinid": np.arange(3, dtype=np.int64),
            "unicode": [
                "\N{LATIN CAPITAL LETTER E}\N{COMBINING CIRCUMFLEX ACCENT}",
                "a \N{GREEK CAPITAL LETTER DELTA} test",
                "クン キン おし.える よ.む くん.ずる",
            ],
            "ascii": ["aa", "bbb", "cccccc"],
            "bytes": [b"aa", b"bbb", b"ccc"],
            "float32": np.array([0.0, 1.0, 2.0], np.float32),
        }
    )

    # We "know" that tiledbsoma will promote string->large_string and
    # binary->large_binary, so generate a table with that latent expectation.
    tbl = pa.Table.from_pandas(df)
    schema = tbl.schema
    for i in range(len(tbl.schema)):
        fld = schema.field(i)
        if fld.type == pa.string():
            schema = schema.set(i, fld.with_type(pa.large_string()))
        if fld.type == pa.binary():
            schema = schema.set(i, fld.with_type(pa.large_binary()))

    return pa.Table.from_pandas(df, schema)


@pytest.fixture
def sample_soma_dataframe(tmp_path, sample_arrow_table):
    sdf = soma.DataFrame(tmp_path.as_posix())
    sdf.create(sample_arrow_table.schema, index_column_names=["soma_joinid"])
    sdf.write(sample_arrow_table)
    assert sdf.exists()
    return sdf


def test_dataframe_unicode_columns(tmp_path, sample_arrow_table):
    """Verify round-trip of unicode in DataFrame value columns"""
    sdf = soma.DataFrame(tmp_path.as_posix())
    sdf.create(sample_arrow_table.schema, index_column_names=["soma_joinid"])
    sdf.write(sample_arrow_table)

    assert sample_arrow_table.schema == sdf.schema
    assert sdf.read().concat().equals(sample_arrow_table)


def test_dataframe_unicode_value_filter(sample_soma_dataframe):
    """Verify that value_filter works correctly"""

    # filter on ascii
    assert sample_soma_dataframe.read(
        value_filter="ascii in ['aa', 'cccccc']"
    ).concat().to_pydict() == {
        "soma_joinid": [0, 2],
        "unicode": [
            "\N{LATIN CAPITAL LETTER E}\N{COMBINING CIRCUMFLEX ACCENT}",
            "クン キン おし.える よ.む くん.ずる",
        ],
        "ascii": ["aa", "cccccc"],
        "bytes": [b"aa", b"ccc"],
        "float32": [0.0, 2.0],
    }

    # filter on unicode, equality
    assert sample_soma_dataframe.read(
        value_filter="unicode == '\N{LATIN CAPITAL LETTER E}\N{COMBINING CIRCUMFLEX ACCENT}'"
    ).concat().to_pydict() == {
        "soma_joinid": [0],
        "unicode": [
            "\N{LATIN CAPITAL LETTER E}\N{COMBINING CIRCUMFLEX ACCENT}",
        ],
        "ascii": ["aa"],
        "bytes": [b"aa"],
        "float32": [0.0],
    }

    # filter on unicode, ordering
    assert sample_soma_dataframe.read(
        value_filter="unicode > 'a'"
    ).concat().to_pydict() == {
        "soma_joinid": [1, 2],
        "unicode": [
            "a \N{GREEK CAPITAL LETTER DELTA} test",
            "クン キン おし.える よ.む くん.ずる",
        ],
        "ascii": ["bbb", "cccccc"],
        "bytes": [b"bbb", b"ccc"],
        "float32": [1.0, 2.0],
    }


# TODO: Remove the `xfail` annotation when TileDB core supports Unicode
# dimensions, aka SOMA index columns.
#
@pytest.mark.xfail
def test_dataframe_unicode_index(tmp_path, sample_arrow_table):
    """Verify round-trip of unicode in DataFrame index columns"""
    sdf = soma.DataFrame(tmp_path.as_posix())
    sdf.create(sample_arrow_table.schema, index_column_names=["unicode"])
    sdf.write(sample_arrow_table)
    assert sdf.read().concat().equals(sample_arrow_table)
