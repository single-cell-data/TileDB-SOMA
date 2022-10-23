import pandas as pd
import pyarrow as pa
import pytest

import tiledbsoma as soma


def test_soma_dataframe_non_indexed(tmp_path):
    sdf = soma.SOMADataFrame(uri=tmp_path.as_posix())

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
    table = sdf.read_all(ids=[1, 2])
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
    # Read by ids
    table = next(sdf.read(ids=pa.array([1, 3])))
    assert table.num_rows == 2

    # We should be getting back the soma_rowid & soma_joinid column as well
    assert table.num_columns == 5

    assert [e.as_py() for e in list(table["soma_rowid"])] == [1, 3]
    assert sorted([e.as_py() for e in list(table["foo"])]) == [20, 40]
    assert sorted([e.as_py() for e in list(table["bar"])]) == [5.2, 7.4]
    assert sorted([e.as_py() for e in list(table["baz"])]) == ["ball", "dog"]

    # ----------------------------------------------------------------
    # Read by value_filter
    table = sdf.read_all(value_filter="foo == 40 or foo == 20")
    assert table.num_rows == 2

    # We should be getting back the soma_rowid & soma_joinid column as well
    assert table.num_columns == 5

    assert [e.as_py() for e in list(table["soma_rowid"])] == [1, 3]
    assert sorted([e.as_py() for e in list(table["foo"])]) == [20, 40]
    assert sorted([e.as_py() for e in list(table["bar"])]) == [5.2, 7.4]
    assert sorted([e.as_py() for e in list(table["baz"])]) == ["ball", "dog"]

    # ----------------------------------------------------------------
    # Read by value_filter
    table = sdf.read_all(value_filter='baz == "ball" or baz == "dog"')
    assert table.num_rows == 2

    # We should be getting back the soma_rowid & soma_joind column as well
    assert table.num_columns == 5

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
            ("soma_rowid", pa.int64()),
            ("soma_joinid", pa.int64()),
            ("A", pa.int64()),
            ("B", pa.float64()),
            ("C", pa.large_string()),
        ]
    )
    sdf = soma.SOMADataFrame(uri=tmp_path.as_posix()).create(schema)
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

    def _check_tbl(tbl, col_names, ids, *, demote):
        assert tbl.num_columns == (
            len(schema.names) if col_names is None else len(col_names)
        )
        assert tbl.num_rows == (n_data if ids is None else len(ids))

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


def test_soma_columns(tmp_path):
    """
    1. soma_joinid/soma_rowid are int64
    2. both will be added by default, if missing in call to create
    3. both are explicit in keys/schema
    4. no other soma_ ids allowed
    """

    A = soma.SOMADataFrame((tmp_path / "A").as_posix())
    A.create(pa.schema([("a", pa.bool_())]))
    assert sorted(A.keys()) == sorted(["a", "soma_rowid", "soma_joinid"])
    assert A.schema.field("soma_joinid").type == pa.int64()
    assert A.schema.field("soma_rowid").type == pa.int64()
    A.delete()

    B = soma.SOMADataFrame((tmp_path / "B").as_posix())
    with pytest.raises(TypeError):
        B.create(pa.schema([("a", pa.bool_()), ("soma_rowid", pa.float32())]))

    C = soma.SOMADataFrame((tmp_path / "C").as_posix())
    with pytest.raises(TypeError):
        C.create(pa.schema([("a", pa.bool_()), ("soma_joinid", pa.float32())]))

    D = soma.SOMADataFrame((tmp_path / "D").as_posix())
    D.create(
        pa.schema(
            [("a", pa.bool_()), ("soma_joinid", pa.int64()), ("soma_rowid", pa.int64())]
        )
    )
    assert sorted(D.keys()) == sorted(["a", "soma_rowid", "soma_joinid"])
    assert D.schema.field("soma_joinid").type == pa.int64()
    assert D.schema.field("soma_rowid").type == pa.int64()
    D.delete()

    E = soma.SOMADataFrame((tmp_path / "E").as_posix())
    with pytest.raises(ValueError):
        E.create(pa.schema([("soma_test", pa.bool_())]))
