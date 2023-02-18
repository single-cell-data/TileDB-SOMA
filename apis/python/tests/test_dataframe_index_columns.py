import numpy as np
import pyarrow as pa
import pytest

import tiledbsoma as soma

# ================================================================
# TODO: SINGLE-INDEX CHECKLIST
# * Index-column names are rows in this table
# * Query types are columns

# .            py-list py-tuple py-slice np-array-untyped np-array-typed pa-array-untyped pa.array-typed
# soma_joinid  .       .        .        .                .              .                .
# string       .       .        .        .                .              .                .
# bytes        .       .        .        .                .              .                .
# int8         .       .        .        .                .              .                .
# int16        .       .        .        .                .              .                .
# int32        .       .        .        .                .              .                .
# int64        .       .        .        .                .              .                .
# uint8        .       .        .        .                .              .                .
# uint16       .       .        .        .                .              .                .
# uint32       .       .        .        .                .              .                .
# uint64       .       .        .        .                .              .                .
# float32      .       .        .        .                .              .                .
# float64      .       .        .        .                .              .                .
# bool         .       .        .        .                .              .                .
# timestamp-s  .       .        .        .                .              .                .
# timestamp-ms .       .        .        .                .              .                .
# timestamp-us .       .        .        .                .              .                .
# timestamp-ns .       .        .        .                .              .                .
# ================================================================

# ================================================================
# TODO: set_dim_ranges as well (slices)
# ================================================================


@pytest.fixture
def arrow_table():
    pydict = {
        "soma_joinid": pa.array([0, 1, 2, 3, 4], pa.int64()),
        "string": ["apple", "ball", "cat", "dog", "egg"],
        "bytes": [b"apple", b"ball", b"cat", b"dog", b"egg"],
        "int8": pa.array([80, 81, 82, 83, 84], pa.int8()),
        "int16": pa.array([1600, 1601, 1602, 1603, 1604], pa.int16()),
        "int32": pa.array([3200, 3201, 3202, 3203, 3204], pa.int32()),
        "int64": pa.array([6400, 6401, 6402, 6403, 6404], pa.int64()),
        "uint8": pa.array([90, 91, 92, 93, 94], pa.uint8()),
        "uint16": pa.array([1610, 1611, 1612, 1613, 1614], pa.uint16()),
        "uint32": pa.array([3210, 3211, 3212, 3213, 3214], pa.uint32()),
        "uint64": pa.array([6410, 6411, 6412, 6413, 6414], pa.uint64()),
        "float32": pa.array([320.5, 321.5, 322.5, 323.5, 324.5], pa.float32()),
        "float64": pa.array([640.5, 641.5, 642.5, 643.5, 644.5], pa.float64()),
        "bool": pa.array([True, True, False, True, False], pa.bool_()),
        "ts": pa.array([1000, 1001, 1002, 1003, 1004], pa.timestamp("s")),
        "tms": pa.array([2000, 2001, 2002, 2003, 2004], pa.timestamp("ms")),
        "tus": pa.array([3000, 3001, 3002, 3003, 3004], pa.timestamp("us")),
        "tns": pa.array([4000, 4001, 4002, 4003, 4004], pa.timestamp("ns")),
    }
    return pa.Table.from_pydict(pydict)


@pytest.mark.parametrize(
    "name,index_column_names,coords,expecteds",
    [
        # Index by soma_joinid
        [
            "soma_joinid/all",
            ["soma_joinid"],
            [],
            {
                "soma_joinid": pa.array([0, 1, 2, 3, 4], pa.int64()),
                "string": pa.array(
                    ["apple", "ball", "cat", "dog", "egg"], pa.large_string()
                ),
            },
        ],
        [
            "soma_joinid/pylist",
            ["soma_joinid"],
            [[0, 2]],
            {
                "soma_joinid": pa.array([0, 2], pa.int64()),
                "string": pa.array(["apple", "cat"], pa.large_string()),
            },
        ],
        [
            "soma_joinid/pytuple",
            ["soma_joinid"],
            [(0, 2)],
            {
                "soma_joinid": pa.array([0, 2], pa.int64()),
                "string": pa.array(["apple", "cat"], pa.large_string()),
            },
        ],
        [
            "soma_joinid/pyslice",
            ["soma_joinid"],
            [slice(0, 2)],
            {
                "soma_joinid": pa.array([0, 1, 2], pa.int64()),
                "string": pa.array(["apple", "ball", "cat"], pa.large_string()),
            },
        ],
        [
            "soma_joinid/np-array-untyped",
            ["soma_joinid"],
            [np.array([0, 2])],
            {
                "soma_joinid": pa.array([0, 2], pa.int64()),
                "string": pa.array(["apple", "cat"], pa.large_string()),
            },
        ],
        [
            "soma_joinid/np-array-typed",
            ["soma_joinid"],
            [np.array([0, 2], np.int64)],
            {
                "soma_joinid": pa.array([0, 2], pa.int64()),
                "string": pa.array(["apple", "cat"], pa.large_string()),
            },
        ],
        [
            "soma_joinid/pa-array-untyped",
            ["soma_joinid"],
            [pa.array([0, 2])],
            {
                "soma_joinid": pa.array([0, 2]),
                "string": pa.array(["apple", "cat"], pa.large_string()),
            },
        ],
        [
            "soma_joinid/pa-array-typed",
            ["soma_joinid"],
            [pa.array([0, 2])],
            {
                "soma_joinid": pa.array([0, 2], pa.int64()),
                "string": pa.array(["apple", "cat"], pa.large_string()),
            },
        ],
        # Index by string
        [
            "string/all",
            ["string"],
            [],
            {
                "soma_joinid": pa.array([0, 1, 2, 3, 4], pa.int64()),
                "string": pa.array(
                    ["apple", "ball", "cat", "dog", "egg"], pa.large_string()
                ),
            },
        ],
        [
            "string/pylist",
            ["string"],
            [["cat", "dog"]],
            {
                "soma_joinid": pa.array([2, 3], pa.int64()),
            },
        ],
        [
            "string/pytuple",
            ["string"],
            [("cat", "dog")],
            {
                "soma_joinid": pa.array([2, 3], pa.int64()),
            },
        ],
        [
            "string/pyslice",
            ["string"],
            [("cat", "dog")],
            {
                "soma_joinid": pa.array([2, 3], pa.int64()),
            },
        ],
        [
            "string/arrow-untyped",
            ["string"],
            [pa.array(["cat", "dog"])],
            {
                "soma_joinid": pa.array([2, 3], pa.int64()),
            },
        ],
        [
            "string/arrow-typed",
            ["string"],
            [pa.array(["cat", "dog"], pa.string())],
            {
                "soma_joinid": pa.array([2, 3], pa.int64()),
            },
        ],
        [
            "string/numpy-untyped",
            ["string"],
            [np.asarray(["cat", "dog"])],
            {
                "soma_joinid": pa.array([2, 3], pa.int64()),
            },
        ],
        [
            "string/numpy-typed",
            ["string"],
            [np.asarray(["cat", "dog"], str)],
            {
                "soma_joinid": pa.array([2, 3], pa.int64()),
            },
        ],
        # Index by int64
        [
            "int64/all",
            ["int64"],
            [],
            {
                "soma_joinid": pa.array([0, 1, 2, 3, 4], pa.int64()),
                "string": pa.array(
                    ["apple", "ball", "cat", "dog", "egg"], pa.large_string()
                ),
            },
        ],
        [
            "int64/arrow",
            ["int64"],
            [pa.array([6402, 6403])],
            {
                "soma_joinid": pa.array([2, 3], pa.int64()),
                "string": pa.array(["cat", "dog"], pa.large_string()),
            },
        ],
        [
            "int64/numpy",
            ["int64"],
            [np.asarray([6402, 6403])],
            {
                "soma_joinid": pa.array([2, 3], pa.int64()),
                "string": pa.array(["cat", "dog"], pa.large_string()),
            },
        ],
        [
            "int64/pylist",
            ["int64"],
            [[6402, 6403]],
            {
                "soma_joinid": pa.array([2, 3], pa.int64()),
                "string": pa.array(["cat", "dog"], pa.large_string()),
            },
        ],
        [
            "int64/pyslice",
            ["int64"],
            [slice(6402, 6403)],
            {
                "soma_joinid": pa.array([2, 3], pa.int64()),
                "string": pa.array(["cat", "dog"], pa.large_string()),
            },
        ],
        # Index by float32
        [
            "float32/pa-array-untyped",
            ["float32"],
            [pa.array([321.5, 323.5])],
            {
                "soma_joinid": pa.array([2, 3], pa.int64()),
            },
        ],
        [
            "float32/pa-array-typed",
            ["float32"],
            [pa.array([322.5, 323.5], pa.float32())],
            {
                "soma_joinid": pa.array([2, 3], pa.int64()),
            },
        ],
        [
            "float32/numpy-untyped",
            ["float32"],
            [np.asarray([322.5, 323.5])],
            {
                "soma_joinid": pa.array([2, 3], pa.int64()),
            },
        ],
        [
            "float32/numpy-typed",
            ["float32"],
            [np.asarray([322.5, 323.5], np.float32)],
            {
                "soma_joinid": pa.array([2, 3], pa.int64()),
            },
        ],
        [
            "float32/pylist",
            ["float32"],
            [[322.5, 323.5]],
            {
                "soma_joinid": pa.array([2, 3], pa.int64()),
            },
        ],
        [
            "float32/pyslice",
            ["float32"],
            [slice(321.5, 323.5)],
            {
                "soma_joinid": pa.array([1, 2, 3], pa.int64()),
            },
        ],
        # Index by float64
        [
            "float64/pa-array-untyped",
            ["float64"],
            [pa.array([641.5, 643.5])],
            {
                "soma_joinid": pa.array([1, 3], pa.int64()),
            },
        ],
        [
            "float64/pa-array-typed",
            ["float64"],
            [pa.array([641.5, 643.5], pa.float64())],
            {
                "soma_joinid": pa.array([1, 3], pa.int64()),
            },
        ],
        [
            "float64/numpy-untyped",
            ["float64"],
            [np.asarray([641.5, 643.5])],
            {
                "soma_joinid": pa.array([1, 3], pa.int64()),
            },
        ],
        [
            "float64/numpy-typed",
            ["float64"],
            [np.asarray([641.5, 643.5], np.float64)],
            {
                "soma_joinid": pa.array([1, 3], pa.int64()),
            },
        ],
        [
            "float64/pylist",
            ["float64"],
            [[641.5, 643.5]],
            {
                "soma_joinid": pa.array([1, 3], pa.int64()),
            },
        ],
        [
            "float64/pyslice",
            ["float64"],
            [slice(641.5, 643.5)],
            {
                "soma_joinid": pa.array([1, 2, 3], pa.int64()),
            },
        ],
        # Index by int64+string
        [
            # TODO: fix flagged bug in _dataframe.py's write method
            "int64+string/arrow",
            ["int64", "string"],
            [pa.array([6402, 6403]), pa.array(["cat", "dog"])],
            {
                "soma_joinid": pa.array([2, 3], pa.int64()),
                "string": pa.array(["cat", "dog"], pa.large_string()),
            },
        ],
        [
            "string+int64/arrow",
            ["string", "int64"],
            [pa.array(["cat", "dog"]), pa.array([6402, 6403])],
            {
                "soma_joinid": pa.array([2, 3], pa.int64()),
                "string": pa.array(["cat", "dog"], pa.large_string()),
            },
        ],
        [
            "string+int64/numpy",
            ["string", "int64"],
            [np.asarray(["cat", "dog"]), np.asarray([6402, 6403])],
            {
                "soma_joinid": pa.array([2, 3], pa.int64()),
                "string": pa.array(["cat", "dog"], pa.large_string()),
            },
        ],
        [
            "string+int64/pylist",
            ["string", "int64"],
            [["cat", "dog"], [6402, 6403]],
            {
                "soma_joinid": pa.array([2, 3], pa.int64()),
                "string": pa.array(["cat", "dog"], pa.large_string()),
            },
        ],
        [
            "string+int64/pytuple",
            ["string", "int64"],
            [("cat", "dog"), (6402, 6403)],
            {
                "soma_joinid": pa.array([2, 3], pa.int64()),
                "string": pa.array(["cat", "dog"], pa.large_string()),
            },
        ],
        # Index by int64+float64+string
        # TODO
    ],
)
def test_types(tmp_path, arrow_table, name, index_column_names, coords, expecteds):
    uri = tmp_path.as_posix()

    soma.DataFrame.create(
        uri,
        schema=arrow_table.schema,
        index_column_names=index_column_names,
    )
    with soma.DataFrame.open(uri, "w") as sdf:
        sdf.write(arrow_table)

    with soma.DataFrame.open(uri, "r") as sdf:
        actual_table = sdf.read(coords=coords).concat()
        for query_column_name, expected_array in expecteds.items():
            actual_array = actual_table[query_column_name].combine_chunks()

            # The output from the first one is easier to read when it fails
            assert actual_array.to_pylist() == expected_array.to_pylist()
            assert actual_array == expected_array
