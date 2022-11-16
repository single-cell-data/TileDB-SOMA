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
            "float32": np.array([0.0, 1.1, 2.2], np.float32),
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


# TODO: Remove the `xfail` annotation when issue TileDB-SOMA#415 is fixed
#
@pytest.mark.xfail
def test_dataframe_unicode(tmp_path, sample_arrow_table):
    """Verify round-trip of unicode in DataFrame attributes"""
    sdf = soma.DataFrame(tmp_path.as_posix())
    sdf.create(sample_arrow_table.schema)
    sdf.write(sample_arrow_table)
    assert sdf.read_all().equals(sample_arrow_table)


# TODO: Remove the `xfail` annotation when issue TileDB-SOMA#415 is fixed
#
@pytest.mark.xfail
def test_indexed_dataframe_unicode_attr(tmp_path, sample_arrow_table):
    """Verify round-trip of unicode in DataFrame value columns"""
    sdf = soma.DataFrame(tmp_path.as_posix())
    sdf.create(sample_arrow_table.schema, index_column_names=["soma_joinid"])
    sdf.write(sample_arrow_table)
    assert sdf.read_all().equals(sample_arrow_table)


# TODO: Remove the `xfail` annotation when issues TileDB-SOMA#415 and TileDB-SOMA#418 are fixed
#
@pytest.mark.xfail
def test_indexed_dataframe_unicode_index(tmp_path, sample_arrow_table):
    """Verify round-trip of unicode in DataFrame index columns"""
    sdf = soma.DataFrame(tmp_path.as_posix())
    sdf.create(sample_arrow_table.schema, index_column_names=["unicode"])
    sdf.write(sample_arrow_table)
    assert sdf.read_all().equals(sample_arrow_table)
