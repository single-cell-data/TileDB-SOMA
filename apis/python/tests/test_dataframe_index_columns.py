import numpy as np
import pandas as pd
import pyarrow as pa
import pytest

import tiledbsoma as soma


@pytest.fixture
def arrow_table():
    pandas_df = pd.DataFrame(
        data={
            "soma_joinid": np.array([0, 1, 2, 3, 4], np.int64),
            # XXX bool
            "string": ["apple", "ball", "cat", "dog", "egg"],
            "bytes": [b"apple", b"ball", b"cat", b"dog", b"egg"],
            "int8": np.array([80, 81, 82, 83, 84], np.int8),
            "int16": np.array([1600, 1601, 1602, 1603, 1604], np.int16),
            "int32": np.array([3200, 3201, 3202, 3203, 3204], np.int32),
            "int64": np.array([6400, 6401, 6402, 6403, 6404], np.int64),
            "uint8": np.array([93, 94, 95, 96, 97], np.uint8),
            "uint16": np.array([1610, 1611, 1612, 1613, 1614], np.uint16),
            "uint32": np.array([3210, 3211, 3212, 3213, 3214], np.uint32),
            "uint64": np.array([6410, 6411, 6412, 6413, 6414], np.uint64),
            "float32": np.array([32.0, 32.1, 32.2, 32.3, 32.4], np.float32),
            "float64": np.array([64.0, 64.1, 64.2, 64.3, 64.4], np.float64),
            # XXX timestamps x 4
        }
    )
    return pa.Table.from_pandas(pandas_df)


@pytest.mark.parametrize(
    "params",
    [
        {
            "index_column_names": ["soma_joinid"],
            "queries": [
                {
                    "coords": [],
                    "expecteds": {
                        "soma_joinid": pa.array([0, 1, 2, 3, 4], pa.int64()),
                        "string": pa.array(
                            ["apple", "ball", "cat", "dog", "egg"], pa.large_string()
                        ),
                    },
                },
                {
                    "coords": [pa.array([0, 2])],
                    "expecteds": {
                        "soma_joinid": pa.array([0, 2], pa.int64()),
                        "string": pa.array(["apple", "cat"], pa.large_string()),
                    },
                },
            ],
        },
        {
            "index_column_names": ["string"],
            "queries": [
                {
                    "coords": [],
                    "expecteds": {
                        "soma_joinid": pa.array([0, 1, 2, 3, 4], pa.int64()),
                        "string": pa.array(
                            ["apple", "ball", "cat", "dog", "egg"], pa.large_string()
                        ),
                    },
                },
                {
                    "coords": [["cat", "dog"]],
                    "expecteds": {
                        "soma_joinid": pa.array([2, 3], pa.int64()),
                    },
                },
                {
                    "coords": [pa.array(["cat", "dog"])],
                    "expecteds": {
                        "soma_joinid": pa.array([2, 3], pa.int64()),
                    },
                },
            ],
        },
        {
            "index_column_names": ["int64"],
            "queries": [
                {
                    "coords": [],
                    "expecteds": {
                        "soma_joinid": pa.array([0, 1, 2, 3, 4], pa.int64()),
                        "string": pa.array(
                            ["apple", "ball", "cat", "dog", "egg"], pa.large_string()
                        ),
                    },
                },
                {
                    "coords": [pa.array([6402, 6403])],
                    "expecteds": {
                        "soma_joinid": pa.array([2, 3], pa.int64()),
                        "string": pa.array(["cat", "dog"], pa.large_string()),
                    },
                },
            ],
        },
        # TODO: many more
        # coordses = {
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
        # }
        #        {
        #            "index_column_names": ["float64"],
        #            "queries": [
        #                {
        #                    "coords": [pa.array([64.1, 64.4])],
        #                    "expecteds": {
        #                        "soma_joinid": pa.array([1, 4], pa.int64()),
        #                    },
        #                },
        #            ],
        #        },
    ],
)
def test_dataframe_index_column_types(tmp_path, arrow_table, params):
    uri = tmp_path.as_posix()

    soma.DataFrame.create(
        uri,
        schema=arrow_table.schema,
        index_column_names=params["index_column_names"],
    )
    with soma.DataFrame.open(uri, "w") as sdf:
        sdf.write(arrow_table)

    for query in params["queries"]:
        coords = query["coords"]
        expecteds = query["expecteds"]
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
