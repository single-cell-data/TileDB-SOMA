import pandas as pd
import pyarrow as pa
import pytest
import tiledb

import tiledbsoma as soma


def test_dataframe_non_indexed(tmp_path):
    sdf = soma.DataFrame(uri=tmp_path.as_posix())

    # Create
    asch = pa.schema(
        [
            ("foo", pa.int32()),
            ("bar", pa.float64()),
            ("baz", pa.large_string()),
        ]
    )
    sdf.create(schema=asch)

    assert sorted(sdf.schema.names) == sorted(
        ["foo", "bar", "baz", "soma_rowid", "soma_joinid"]
    )
    assert sorted(sdf.keys()) == sorted(sdf.schema.names)

    # ----------------------------------------------------------------
    # Write
    for _i in range(3):
        pydict = {}
        pydict["soma_rowid"] = [0, 1, 2, 3, 4]
        pydict["soma_joinid"] = [101, 102, 103, 104, 105]
        pydict["foo"] = [10, 20, 30, 40, 50]
        pydict["bar"] = [4.1, 5.2, 6.3, 7.4, 8.5]
        pydict["baz"] = ["apple", "ball", "cat", "dog", "egg"]
        rb = pa.Table.from_pydict(pydict)
        sdf.write(rb)

    # ----------------------------------------------------------------
    # Read all
    table = sdf.read_all()
    assert table.num_rows == 5

    # We should be getting back the soma_rowid & soma_joinid column as well
    assert table.num_columns == 5

    assert [e.as_py() for e in list(table["soma_rowid"])] == [0, 1, 2, 3, 4]
    assert [e.as_py() for e in list(table["soma_joinid"])] == [101, 102, 103, 104, 105]
    assert [e.as_py() for e in list(table["foo"])] == pydict["foo"]
    assert [e.as_py() for e in list(table["bar"])] == pydict["bar"]
    assert [e.as_py() for e in list(table["baz"])] == pydict["baz"]

    # ----------------------------------------------------------------
    # Read by ids
    table = sdf.read_all(ids=slice(1, 2))
    assert table.num_rows == 2

    # We should be getting back the soma_rowid column as well
    # If sparse dataframe:
    assert table.num_columns == 5
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

    # We should be getting back the soma_rowid & soma_joinid column as well
    assert table.num_columns == 5

    assert [e.as_py() for e in list(table["soma_rowid"])] == [1, 2]
    assert sorted([e.as_py() for e in list(table["foo"])]) == [20, 30]
    assert sorted([e.as_py() for e in list(table["bar"])]) == [5.2, 6.3]
    assert sorted([e.as_py() for e in list(table["baz"])]) == ["ball", "cat"]

    # ----------------------------------------------------------------
    # Read by value_filter
    table = sdf.read_all(value_filter="foo == 40 or foo == 20")
    assert table.num_rows == 5

    # We should be getting back the soma_rowid & soma_joinid column as well
    assert table.num_columns == 5

    with sdf._tiledb_open() as A:
        mask = A.attr("soma_joinid").fill

    retain_flags = table["soma_joinid"] != mask
    filtered_table = table.filter(retain_flags)

    assert [e.as_py() for e in list(filtered_table["soma_rowid"])] == [1, 3]
    assert sorted([e.as_py() for e in list(filtered_table["foo"])]) == [20, 40]
    assert sorted([e.as_py() for e in list(filtered_table["bar"])]) == [5.2, 7.4]
    assert sorted([e.as_py() for e in list(filtered_table["baz"])]) == ["ball", "dog"]

    # ----------------------------------------------------------------
    # Read by value_filter
    table = sdf.read_all(value_filter='baz == "ball" or baz == "dog"')
    assert table.num_rows == 5

    # We should be getting back the soma_rowid & soma_joind column as well
    assert table.num_columns == 5

    with sdf._tiledb_open() as A:
        mask = A.attr("soma_joinid").fill

    retain_flags = table["soma_joinid"] != mask
    filtered_table = table.filter(retain_flags)

    assert sorted([e.as_py() for e in list(filtered_table["soma_joinid"])]) == [
        102,
        104,
    ]
    assert sorted([e.as_py() for e in list(filtered_table["foo"])]) == [20, 40]
    assert sorted([e.as_py() for e in list(filtered_table["bar"])]) == [5.2, 7.4]
    assert sorted([e.as_py() for e in list(filtered_table["baz"])]) == ["ball", "dog"]


@pytest.fixture
def simple_data_frame(tmp_path):
    """
    A pytest fixture which creates a simple DataFrame for use in tests below.
    """
    schema = pa.schema(
        [
            ("soma_rowid", pa.int64()),
            ("soma_joinid", pa.int64()),
            ("A", pa.int64()),
            ("B", pa.float64()),
            ("C", pa.large_string()),
        ]
    )
    sdf = soma.DataFrame(uri=tmp_path.as_posix()).create(schema)
    data = {
        "soma_rowid": [0, 1, 2, 3],
        "soma_joinid": [100, 200, 300, 400],
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
        0,
        slice(1, 3),
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
def test_DataFrame_read_column_names(simple_data_frame, ids, col_names):
    """
    Issue #312 - `column_names` parameter not correctly handled.

    While the bug report was only against DataFrame.read,this
    test covers all of the read* methods.
    """

    schema, sdf, n_data = simple_data_frame
    assert sdf.exists()

    def _check_tbl(tbl, col_names, ids, *, demote):
        assert tbl.num_columns == (
            len(schema.names) if col_names is None else len(col_names)
        )
        assert tbl.num_rows == (n_data if ids is None else (ids.stop - ids.start + 1))

        if demote:
            assert tbl.schema == pa.schema(
                [
                    pa.field(schema.field(f).name, pa.string())
                    if schema.field(f).type == pa.large_string()
                    else schema.field(f)
                    for f in (col_names if col_names is not None else schema.names)
                ]
            )
        else:
            assert tbl.schema == pa.schema(
                [
                    schema.field(f)
                    for f in (col_names if col_names is not None else schema.names)
                ]
            )

    # TileDB ASCII -> Arrow large_string
    _check_tbl(
        sdf.read_all(ids=ids, column_names=col_names),
        col_names,
        ids,
        demote=False,
    )

    _check_tbl(
        sdf.read_all(column_names=col_names),
        col_names,
        None,
        demote=False,
    )

    # TileDB ASCII -> Pandas string -> Arrow string (not large_string)
    _check_tbl(
        pa.Table.from_pandas(
            pd.concat(sdf.read_as_pandas(ids=ids, column_names=col_names))
        ),
        col_names,
        ids,
        demote=True,
    )

    _check_tbl(
        pa.Table.from_pandas(sdf.read_as_pandas_all(column_names=col_names)),
        col_names,
        None,
        demote=True,
    )


def test_empty_dataframe(tmp_path):
    a = soma.DataFrame((tmp_path / "A").as_posix())
    a.create(pa.schema([("a", pa.bool_())]))
    # Must not throw
    ta1 = next(a.read())
    ta2 = a.read_all()
    tp1 = next(a.read_as_pandas())
    tp2 = a.read_as_pandas_all()

    assert len(ta1) == 1
    assert len(ta2) == 1
    assert len(tp1) == 1
    assert len(tp2) == 1
    assert isinstance(tp2, pd.DataFrame)

    with tiledb.open(a.uri) as A:
        jid_mask = A.attr("soma_joinid").fill
    assert ta1["soma_joinid"][0].as_py() == jid_mask
    assert ta2["soma_joinid"][0].as_py() == jid_mask
    assert tp1["soma_joinid"][0] == jid_mask
    assert tp2["soma_joinid"][0] == jid_mask


def test_columns(tmp_path):
    """
    1. soma_joinid/soma_rowid are int64
    2. both will be added by default, if missing in call to create
    3. both are explicit in keys/schema
    4. no other soma_ ids allowed
    """

    A = soma.DataFrame((tmp_path / "A").as_posix())
    A.create(pa.schema([("a", pa.bool_())]))
    assert sorted(A.keys()) == sorted(["a", "soma_rowid", "soma_joinid"])
    assert A.schema.field("soma_joinid").type == pa.int64()
    assert A.schema.field("soma_rowid").type == pa.int64()
    A.delete()

    B = soma.DataFrame((tmp_path / "B").as_posix())
    with pytest.raises(TypeError):
        B.create(pa.schema([("a", pa.bool_()), ("soma_rowid", pa.float32())]))

    C = soma.DataFrame((tmp_path / "C").as_posix())
    with pytest.raises(TypeError):
        C.create(pa.schema([("a", pa.bool_()), ("soma_joinid", pa.float32())]))

    D = soma.DataFrame((tmp_path / "D").as_posix())
    D.create(
        pa.schema(
            [("a", pa.bool_()), ("soma_joinid", pa.int64()), ("soma_rowid", pa.int64())]
        )
    )
    assert sorted(D.keys()) == sorted(["a", "soma_rowid", "soma_joinid"])
    assert D.schema.field("soma_joinid").type == pa.int64()
    assert D.schema.field("soma_rowid").type == pa.int64()
    D.delete()

    E = soma.DataFrame((tmp_path / "E").as_posix())
    with pytest.raises(ValueError):
        E.create(pa.schema([("soma_test", pa.bool_())]))
