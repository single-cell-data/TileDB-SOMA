import pyarrow as pa
import pytest

import tiledbsoma as soma

# ================================================================
# TODO: SINGLE-INDEX CHECKLIST
# * Index-column names are columns in this table
# * Query types are rows
# * Example: "int32-pa" means:
#   o  Set index_column_names to ["int32"]
#   o  Query using a length-one list of pyarrow array
#
# .               soma_joinid string bytes int8 int16 int32 int64 uint8 uint16 uint32 uint64 float32 float64
# soma_joinid-py  .           .      .     .    .     .     .     .     .      .      .      .       .
# soma_joinid-pa  .           .      .     .    .     .     .     .     .      .      .      .       .
# soma_joinid-np  .           .      .     .    .     .     .     .     .      .      .      .       .
# string-py       .           .      .     .    .     .     .     .     .      .      .      .       .
# string-pas      .           .      .     .    .     .     .     .     .      .      .      .       .
# string-pal      .           .      .     .    .     .     .     .     .      .      .      .       .
# string-np       .           .      .     .    .     .     .     .     .      .      .      .       .
# bytes-py        .           .      .     .    .     .     .     .     .      .      .      .       .
# bytes-pas       .           .      .     .    .     .     .     .     .      .      .      .       .
# bytes-pal       .           .      .     .    .     .     .     .     .      .      .      .       .
# bytes-np        .           .      .     .    .     .     .     .     .      .      .      .       .
# int-py          .           .      .     .    .     .     .     .     .      .      .      .       .
# int8-pa         .           .      .     .    .     .     .     .     .      .      .      .       .
# int8-np         .           .      .     .    .     .     .     .     .      .      .      .       .
# int16-pa        .           .      .     .    .     .     .     .     .      .      .      .       .
# int16-np        .           .      .     .    .     .     .     .     .      .      .      .       .
# int32-pa        .           .      .     .    .     .     .     .     .      .      .      .       .
# int32-np        .           .      .     .    .     .     .     .     .      .      .      .       .
# int64-pa        .           .      .     .    .     .     .     .     .      .      .      .       .
# int64-np        .           .      .     .    .     .     .     .     .      .      .      .       .
# uint8-pa        .           .      .     .    .     .     .     .     .      .      .      .       .
# uint8-np        .           .      .     .    .     .     .     .     .      .      .      .       .
# uint16-pa       .           .      .     .    .     .     .     .     .      .      .      .       .
# uint16-np       .           .      .     .    .     .     .     .     .      .      .      .       .
# uint32-pa       .           .      .     .    .     .     .     .     .      .      .      .       .
# uint32-np       .           .      .     .    .     .     .     .     .      .      .      .       .
# uint64-pa       .           .      .     .    .     .     .     .     .      .      .      .       .
# uint64-np       .           .      .     .    .     .     .     .     .      .      .      .       .
# float-py        .           .      .     .    .     .     .     .     .      .      .      .       .
# float32-pa      .           .      .     .    .     .     .     .     .      .      .      .       .
# float32-np      .           .      .     .    .     .     .     .     .      .      .      .       .
# float64-pa      .           .      .     .    .     .     .     .     .      .      .      .       .
# float64-np      .           .      .     .    .     .     .     .     .      .      .      .       .
# ================================================================

# ================================================================
# TODO: MULTI-INDEX CHECKLIST
# * Index-column names are columns in this table
# * Query types are rows
# * Example: "string-int32-pa" means:
#   o  Set index_column_names to ["int32", "string"]
#   o  Query using a list of pyarrow arrays
#
# .               int32+string string+int32
# int32-string-py .            .
# string-int32-py .            .
# int32-string-pa .            .
# string-int32-pa .            .
# int32-string-np .            .
# string-int32-np .            .
#
# .                       int32+float64+string string+int32+float
# int32-string-float64-py .                    .
# string-int32-float64-py .                    .
# int32-string-float64-pa .                    .
# float64-string-int32-pa .                    .
# int32-float64-string-np .                    .
# string-int32-float64-np .                    .
# ================================================================

# TODO: set_dim_ranges as well (slices)

# TODO: many more
#    "soma_joinid": [[0,2]],
#
#    "stringpy":   [["apple", "egg"]],
#    "stringnp": [np.array(["apple", "egg"])],
#    "stringpa": [pa.array(["apple", "egg"])],
#
#    "bytespy":   [[b"apple", b"egg"]],
#    "bytesnp": [np.array([b"apple", b"egg"])],
#    "bytespa": [pa.array([b"apple", b"egg"])],
#
#    "intpy":   [[7, 6]],
#
#    "int8np":  [np.array([87, 86], np.int8)],
#    "int8pa":  [pa.array([87, 86], pa.int8())],
#    "int16np": [np.array([1607, 1606], np.int16)],
#    "int16pa": [pa.array([1607, 1606], pa.int16())],
#    "int32np": [np.array([3207, 3206], np.int32)],
#    "int32pa": [pa.array([3207, 3206], pa.int32())],
#    "int64np": [np.array([6407, 6406], np.int64)],
#    "int64pa": [pa.array([6407, 6406], pa.int64())],
#
#    "uint8np":  [np.array([97, 96], np.uint8)],
#    "uint8pa":  [pa.array([97, 96], pa.uint8())],
#    "uint16np": [np.array([1617, 1616], np.uint16)],
#    "uint16pa": [pa.array([1617, 1616], pa.uint16())],
#    "uint32np": [np.array([3217, 3216], np.uint32)],
#    "uint32pa": [pa.array([3217, 3216], pa.uint32())],
#    "uint64np": [np.array([6417, 6416], np.uint64)],
#    "uint64pa": [pa.array([6417, 6416], pa.uint64())],
#
#    "floatpy": [[32.3, 64.4]],
#    "float32np": [np.array([32.3, 32.4], np.float32)],
#    "float32pa": [pa.array([32.3, 32.4], pa.float32())],
#    "float64np": [np.array([64.2, 64.3], np.float64)],
#    "float64pa": [pa.array([64.2, 64.3], pa.float64())],


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
        "uint8": pa.array([93, 94, 95, 96, 97], pa.uint8()),
        "uint16": pa.array([1610, 1611, 1612, 1613, 1614], pa.uint16()),
        "uint32": pa.array([3210, 3211, 3212, 3213, 3214], pa.uint32()),
        "uint64": pa.array([6410, 6411, 6412, 6413, 6414], pa.uint64()),
        "float32": pa.array([32.0, 32.1, 32.2, 32.3, 32.4], pa.float32()),
        "float64": pa.array([64.0, 64.1, 64.2, 64.3, 64.4], pa.float64()),
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
            "soma_joinid/arrow",
            ["soma_joinid"],
            [pa.array([0, 2])],
            {
                "soma_joinid": pa.array([0, 2], pa.int64()),
                "string": pa.array(["apple", "cat"], pa.large_string()),
            },
        ],
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
            "string/arrow",
            ["string"],
            [pa.array(["cat", "dog"])],
            {
                "soma_joinid": pa.array([2, 3], pa.int64()),
            },
        ],
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
        # [
        # TODO: fix flagged bug in _dataframe.py's write method
        #        {
        #            "index_column_names": ["int64", "string"],
        #            "queries": [
        #                {
        #                    "coords": [pa.array([6402, 6403]), pa.array(["cat", "dog"])],
        #                    "expecteds": {
        #                        "soma_joinid": pa.array([2, 3], pa.int64()),
        #                        "string": pa.array(["cat", "dog"], pa.large_string()),
        #                    },
        #                },
        #            ],
        #        },
        # ],
        [
            "string+int64/arrow",
            ["string", "int64"],
            [pa.array(["cat", "dog"]), pa.array([6402, 6403])],
            {
                "soma_joinid": pa.array([2, 3], pa.int64()),
                "string": pa.array(["cat", "dog"], pa.large_string()),
            },
        ],
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

            # print()
            # print("LEFT")
            # print(type(actual_array))
            # print(actual_array)
            # print("RIGHT")
            # print(type(expected_array))
            # print(expected_array)

            # The output from the first one is easier to read when it fails
            assert actual_array.to_pylist() == expected_array.to_pylist()
            assert actual_array == expected_array
