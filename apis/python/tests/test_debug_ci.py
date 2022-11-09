from typing import List

import pyarrow as pa
import pytest

import tiledbsoma as soma


def make_multiply_indexed_dataframe_debug_ci(tmp_path, index_column_names: List[str]):
    """
    Creates a variably-indexed IndexedDataFrame for use in tests below.
    """
    schema = pa.schema(
        [
            # TO DO: Support non-int index types when we have non-int index support
            # in libtiledbsoma's SOMAReader. See also
            # https://github.com/single-cell-data/TileDB-SOMA/issues/418
            # https://github.com/single-cell-data/TileDB-SOMA/issues/419
            ("index1", pa.int64()),
            ("index2", pa.int64()),
            ("index3", pa.int64()),
            ("index4", pa.int64()),
            ("soma_joinid", pa.int64()),
            ("A", pa.int64()),
        ]
    )

    sidf = soma.IndexedDataFrame(uri=tmp_path.as_posix())
    sidf.create(schema=schema, index_column_names=index_column_names)

    data = {
        "index1": [0, 1, 2, 3, 4, 5],
        "index2": [400, 400, 500, 500, 600, 600],
        "index3": [0, 1, 0, 1, 0, 1],
        "index4": [1000, 2000, 1000, 1000, 1000, 1000],
        "soma_joinid": [10, 11, 12, 13, 14, 15],
        "A": [10, 11, 12, 13, 14, 15],
    }

    n_data = len(data["index1"])
    rb = pa.Table.from_pydict(data)
    sidf.write(rb)
    return (schema, sidf, n_data)


@pytest.mark.parametrize(
    "io",
    [
        # 2D: indexing list is None
        {
            "index_column_names": ["index2", "index3"],
            "ids": None,
            "A": [10, 11, 12, 13, 14, 15],
            "throws": None,
        },
        # 2D: indexing slot is None
        {
            "index_column_names": ["index2", "index3"],
            "ids": [None, None],
            "A": [10, 11, 12, 13, 14, 15],
            "throws": None,
        },
    ],
)
def test_debug_ci(tmp_path, io):
    """Test various ways of indexing on read"""

    schema, sidf, n_data = make_multiply_indexed_dataframe_debug_ci(
        tmp_path, io["index_column_names"]
    )
    assert sidf.exists()

    col_names = ["A"]

    if io["throws"] is not None:
        with pytest.raises(io["throws"]):
            next(sidf.read(ids=io["ids"], column_names=col_names))
    else:
        table = next(sidf.read(ids=io["ids"], column_names=col_names))
        assert table["A"].to_pylist() == io["A"]

    if io["throws"] is not None:
        with pytest.raises(io["throws"]):
            next(sidf.read_as_pandas(ids=io["ids"], column_names=col_names))
    else:
        table = next(sidf.read_as_pandas(ids=io["ids"], column_names=col_names))
        assert table["A"].to_list() == io["A"]

    sidf.delete()
