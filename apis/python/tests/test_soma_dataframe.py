import pandas as pd
import pyarrow as pa
import pytest

import tiledbsoma as t


def test_soma_dataframe_non_indexed(tmp_path):
    sdf = t.SOMADataFrame(uri=tmp_path.as_posix())

    # Create
    asch = pa.schema(
        [
            ("foo", pa.int32()),
            ("bar", pa.float64()),
            ("baz", pa.large_string()),
        ]
    )
    sdf.create(schema=asch)

    # ----------------------------------------------------------------
    # Write
    for _i in range(3):
        pydict = {}
        pydict["soma_rowid"] = [0, 1, 2, 3, 4]
        pydict["foo"] = [10, 20, 30, 40, 50]
        pydict["bar"] = [4.1, 5.2, 6.3, 7.4, 8.5]
        pydict["baz"] = ["apple", "ball", "cat", "dog", "egg"]
        rb = pa.Table.from_pydict(pydict)
        sdf.write(rb)

    # ----------------------------------------------------------------
    # Read all
    table = sdf.read_all()
    assert table.num_rows == 5

    # We should be getting back the soma_rowid column as well
    # If sparse dataframe:
    assert table.num_columns == 4
    # If dense dataframe:
    # assert table.num_columns == 3

    # TODO assert [e.as_py() for e in list(table['soma_rowid'])] == [0,1,2,3,4]
    assert [e.as_py() for e in list(table["foo"])] == pydict["foo"]
    assert [e.as_py() for e in list(table["bar"])] == pydict["bar"]
    assert [e.as_py() for e in list(table["baz"])] == pydict["baz"]

    # ----------------------------------------------------------------
    # Read by ids
    table = sdf.read_all(ids=[1, 2])
    assert table.num_rows == 2

    # We should be getting back the soma_rowid column as well
    # If sparse dataframe:
    assert table.num_columns == 4
    # If dense dataframe:
    # assert table.num_columns == 3

    # TODO assert [e.as_py() for e in list(table['soma_rowid'])] == [0,1,2,3,4]
    assert sorted([e.as_py() for e in list(table["foo"])]) == [20, 30]
    assert sorted([e.as_py() for e in list(table["bar"])]) == [5.2, 6.3]
    assert sorted([e.as_py() for e in list(table["baz"])]) == ["ball", "cat"]

    # ----------------------------------------------------------------
    # Read by ids
    table = sdf.read_all(ids=slice(1, 2))
    assert table.num_rows == 2

    # We should be getting back the soma_rowid column as well
    # If sparse dataframe:
    assert table.num_columns == 4
    # If dense dataframe:
    # assert table.num_columns == 3

    # TODO assert [e.as_py() for e in list(table['soma_rowid'])] == [0,1,2,3,4]
    assert sorted([e.as_py() for e in list(table["foo"])]) == [20, 30]
    assert sorted([e.as_py() for e in list(table["bar"])]) == [5.2, 6.3]
    assert sorted([e.as_py() for e in list(table["baz"])]) == ["ball", "cat"]

    # ----------------------------------------------------------------
    # Read by value_filter
    table = sdf.read_all(value_filter="foo == 40 or foo == 20")
    assert table.num_rows == 2

    # We should be getting back the soma_rowid column as well
    # If sparse dataframe:
    assert table.num_columns == 4
    # If dense dataframe:
    # assert table.num_columns == 3

    # TODO assert [e.as_py() for e in list(table['soma_rowid'])] == [0,1,2,3,4]
    assert sorted([e.as_py() for e in list(table["foo"])]) == [20, 40]
    assert sorted([e.as_py() for e in list(table["bar"])]) == [5.2, 7.4]
    assert sorted([e.as_py() for e in list(table["baz"])]) == ["ball", "dog"]

    # ----------------------------------------------------------------
    # Read by value_filter
    table = sdf.read_all(value_filter='baz == "ball" or baz == "dog"')
    assert table.num_rows == 2

    # We should be getting back the soma_rowid column as well
    # If sparse dataframe:
    assert table.num_columns == 4
    # If dense dataframe:
    # assert table.num_columns == 3

    # TODO assert [e.as_py() for e in list(table['soma_rowid'])] == [0,1,2,3,4]
    assert sorted([e.as_py() for e in list(table["foo"])]) == [20, 40]
    assert sorted([e.as_py() for e in list(table["bar"])]) == [5.2, 7.4]
    assert sorted([e.as_py() for e in list(table["baz"])]) == ["ball", "dog"]


@pytest.fixture
def simple_soma_data_frame(tmp_path):
    """
    A pytest fixture which creates a simple SOMADataFrame for use in tests below.
    """
    schema = pa.schema(
        [
            ("soma_rowid", pa.uint64()),
            ("A", pa.int64()),
            ("B", pa.float64()),
            ("C", pa.large_string()),
        ]
    )
    sdf = t.SOMADataFrame(uri=tmp_path.as_posix())

    # TODO: see issue #324 - create will not accept 'soma_rowid' in the schema,
    # even if the type of that field is correctly specified as a uint64.
    sdf.create(schema=schema.remove(schema.get_field_index("soma_rowid")))

    data = {
        "soma_rowid": [0, 1, 2, 3],
        "A": [10, 11, 12, 13],
        "B": [100.1, 200.2, 300.3, 400.4],
        "C": ["this", "is", "a", "test"],
    }
    n_data = len(data["soma_rowid"])
    rb = pa.Table.from_pydict(data)
    sdf.write(rb)
    yield (schema, sdf, n_data)
    sdf.delete()


@pytest.mark.parametrize(
    "ids",
    [
        None,
        [
            0,
        ],
        [1, 3],
    ],
)
@pytest.mark.parametrize(
    "col_names",
    [
        ["A"],
        ["B"],
        ["A", "B"],
        ["soma_rowid"],
        ["soma_rowid", "A", "B", "C"],
        None,
    ],
)
def test_SOMADataFrame_read_column_names(simple_soma_data_frame, ids, col_names):
    """
    Issue #312 - `column_names` parameter not correctly handled.

    While the bug report was only against SOMADataFrame.read,this
    test covers all of the read* methods.
    """

    schema, sdf, n_data = simple_soma_data_frame
    assert sdf.exists()

    def _check_tbl(tbl, col_names, ids):
        assert tbl.num_columns == (
            len(schema.names) if col_names is None else len(col_names)
        )
        assert tbl.num_rows == (n_data if ids is None else len(ids))
        assert tbl.schema == pa.schema(
            [
                schema.field(f)
                for f in (col_names if col_names is not None else schema.names)
            ]
        )

    _check_tbl(
        sdf.read_all(ids=ids, column_names=col_names),
        col_names,
        ids,
    )
    _check_tbl(
        sdf.read_all(column_names=col_names),
        col_names,
        None,
    )
    _check_tbl(
        pa.Table.from_pandas(
            pd.concat(sdf.read_as_pandas(ids=ids, column_names=col_names))
        ),
        col_names,
        ids,
    )
    _check_tbl(
        pa.Table.from_pandas(sdf.read_as_pandas_all(column_names=col_names)),
        col_names,
        None,
    )
