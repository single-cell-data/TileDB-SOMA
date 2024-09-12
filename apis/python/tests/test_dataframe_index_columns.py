import numpy as np
import pyarrow as pa
import pytest

import tiledbsoma as soma


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
        "tss": pa.array(
            [946684800, 946684801, 946684802, 946684803, 946684804], pa.timestamp("s")
        ),
        "tsms": pa.array(
            [946684800000, 946684800001, 946684800002, 946684800003, 946684800004],
            pa.timestamp("ms"),
        ),
        "tsus": pa.array(
            [
                946684800000000,
                946684800000001,
                946684800000002,
                946684800000003,
                946684800000004,
            ],
            pa.timestamp("us"),
        ),
        "tsns": pa.array(
            [
                946684800000000000,
                946684800000000001,
                946684800000000002,
                946684800000000003,
                946684800000000004,
            ],
            pa.timestamp("ns"),
        ),
    }
    return pa.Table.from_pydict(pydict)


@pytest.mark.parametrize(
    "name,index_column_names,domain,coords,expecteds",
    [
        # Index by soma_joinid
        [
            "SOMA_JOINID-ALL",
            ["soma_joinid"],
            None,
            [],
            "default01234",
        ],
        [
            "soma_joinid-all-shaped",
            ["soma_joinid"],
            [[0, 10]],
            [],
            "default01234",
        ],
        [
            "soma_joinid-py-list",
            ["soma_joinid"],
            None,
            [[0, 2]],
            {
                "soma_joinid": pa.array([0, 2], pa.int64()),
                "string": pa.array(["apple", "cat"], pa.large_string()),
            },
        ],
        [
            "soma_joinid-py-tuple",
            ["soma_joinid"],
            None,
            [(0, 2)],
            {
                "soma_joinid": pa.array([0, 2], pa.int64()),
                "string": pa.array(["apple", "cat"], pa.large_string()),
            },
        ],
        [
            "soma_joinid-py-slice",
            ["soma_joinid"],
            None,
            [slice(0, 2)],
            {
                "soma_joinid": pa.array([0, 1, 2], pa.int64()),
                "string": pa.array(["apple", "ball", "cat"], pa.large_string()),
            },
        ],
        [
            "soma_joinid-py-left-none-slice",
            ["soma_joinid"],
            None,
            [slice(None, 2)],
            {
                "soma_joinid": pa.array([0, 1, 2], pa.int64()),
                "string": pa.array(["apple", "ball", "cat"], pa.large_string()),
            },
        ],
        [
            "soma_joinid-py-right-none-slice",
            ["soma_joinid"],
            None,
            [slice(2, None)],
            {
                "soma_joinid": pa.array([2, 3, 4], pa.int64()),
                "string": pa.array(["cat", "dog", "egg"], pa.large_string()),
            },
        ],
        [
            "soma_joinid-py-both-none-slice",
            ["soma_joinid"],
            None,
            [slice(None, None)],
            "default01234",
        ],
        [
            "soma_joinid-np-array-untyped",
            ["soma_joinid"],
            None,
            [np.array([0, 2])],
            {
                "soma_joinid": pa.array([0, 2], pa.int64()),
                "string": pa.array(["apple", "cat"], pa.large_string()),
            },
        ],
        [
            "soma_joinid-np-array-typed",
            ["soma_joinid"],
            None,
            [np.array([0, 2], np.int64)],
            {
                "soma_joinid": pa.array([0, 2], pa.int64()),
                "string": pa.array(["apple", "cat"], pa.large_string()),
            },
        ],
        [
            "soma_joinid-pa-array-untyped",
            ["soma_joinid"],
            None,
            [pa.array([0, 2])],
            {
                "soma_joinid": pa.array([0, 2]),
                "string": pa.array(["apple", "cat"], pa.large_string()),
            },
        ],
        [
            "soma_joinid-pa-array-typed",
            ["soma_joinid"],
            None,
            [pa.array([0, 2])],
            {
                "soma_joinid": pa.array([0, 2], pa.int64()),
                "string": pa.array(["apple", "cat"], pa.large_string()),
            },
        ],
        # Index by string
        [
            "STRING-ALL",
            ["string"],
            None,
            [],
            "default01234",
        ],
        [
            "string-py-list",
            ["string"],
            None,
            [["cat", "dog"]],
            "default23",
        ],
        [
            "string-py-tuple",
            ["string"],
            None,
            [("cat", "dog")],
            "default23",
        ],
        [
            "string-py-slice",
            ["string"],
            None,
            [("cat", "dog")],
            "default23",
        ],
        [
            "string-pa-array-untyped",
            ["string"],
            None,
            [pa.array(["cat", "dog"])],
            "default23",
        ],
        [
            "string-pa-array-typed",
            ["string"],
            None,
            [pa.array(["cat", "dog"], pa.string())],
            "default23",
        ],
        [
            "string-np-array-untyped",
            ["string"],
            None,
            [np.asarray(["cat", "dog"])],
            "default23",
        ],
        [
            "string-np-array-typed",
            ["string"],
            None,
            [np.asarray(["cat", "dog"], str)],
            "default23",
        ],
        # Index by bytes
        [
            "BYTES-ALL",
            ["bytes"],
            None,
            [],
            "default01234",
        ],
        [
            "bytes-py-list",
            ["bytes"],
            None,
            [[b"cat", b"dog"]],
            "default23",
        ],
        [
            "bytes-py-tuple",
            ["bytes"],
            None,
            [(b"cat", b"dog")],
            "default23",
        ],
        [
            "bytes-py-slice",
            ["bytes"],
            None,
            [(b"cat", b"dog")],
            "default23",
        ],
        [
            "bytes-pa-array-untyped",
            ["bytes"],
            None,
            [pa.array([b"cat", b"dog"])],
            "default23",
        ],
        [
            "bytes-pa-array-typed",
            ["bytes"],
            None,
            [pa.array([b"cat", b"dog"], pa.binary())],
            "default23",
        ],
        [
            "bytes-np-array-untyped",
            ["bytes"],
            None,
            [np.asarray([b"cat", b"dog"])],
            "default23",
        ],
        [
            "bytes-np-array-typed",
            ["bytes"],
            None,
            [np.asarray([b"cat", b"dog"], bytes)],
            "default23",
        ],
        # Index by int64
        [
            "INT64-ALL",
            ["int64"],
            None,
            [],
            "default01234",
        ],
        [
            "int64-all-shaped",
            ["int64"],
            [[6400, 6500]],
            [],
            "default01234",
        ],
        [
            "int64-py-list",
            ["int64"],
            None,
            [[6402, 6403]],
            "default23",
        ],
        [
            "int64-py-tuple",
            ["int64"],
            None,
            [[6402, 6403]],
            "default23",
        ],
        [
            "int64-py-slice",
            ["int64"],
            None,
            [slice(6402, 6403)],
            "default23",
        ],
        [
            "int64-py-left-none-slice",
            ["int64"],
            None,
            [slice(None, 6402)],
            {
                "soma_joinid": pa.array([0, 1, 2], pa.int64()),
                "string": pa.array(["apple", "ball", "cat"], pa.large_string()),
            },
        ],
        [
            "int64-py-right-none-slice",
            ["int64"],
            None,
            [slice(6402, None)],
            {
                "soma_joinid": pa.array([2, 3, 4], pa.int64()),
                "string": pa.array(["cat", "dog", "egg"], pa.large_string()),
            },
        ],
        [
            "int64-py-both-none-slice",
            ["int64"],
            None,
            [slice(None, None)],
            "default01234",
        ],
        [
            "int64-numpy-untyped",
            ["int64"],
            None,
            [np.asarray([6402, 6403])],
            "default23",
        ],
        [
            "int64-numpy-typed",
            ["int64"],
            None,
            [np.asarray([6402, 6403], np.int64)],
            "default23",
        ],
        [
            "int64-pa-array-untyped",
            ["int64"],
            None,
            [pa.array([6402, 6403])],
            "default23",
        ],
        [
            "int64-pa-array-typed",
            ["int64"],
            None,
            [pa.array([6402, 6403], pa.int64())],
            "default23",
        ],
        # Index by int32
        [
            "INT32-ALL",
            ["int32"],
            None,
            [],
            "default01234",
        ],
        [
            "int32-all-shaped",
            ["int32"],
            [[3200, 3300]],
            [],
            "default01234",
        ],
        [
            "int32-py-list",
            ["int32"],
            None,
            [[3202, 3203]],
            "default23",
        ],
        [
            "int32-py-tuple",
            ["int32"],
            None,
            [[3202, 3203]],
            "default23",
        ],
        [
            "int32-py-slice",
            ["int32"],
            None,
            [slice(3202, 3203)],
            "default23",
        ],
        [
            "int32-py-left-none-slice",
            ["int32"],
            None,
            [slice(None, 3202)],
            {
                "soma_joinid": pa.array([0, 1, 2], pa.int64()),
                "string": pa.array(["apple", "ball", "cat"], pa.large_string()),
            },
        ],
        [
            "int32-py-right-none-slice",
            ["int32"],
            None,
            [slice(3202, None)],
            {
                "soma_joinid": pa.array([2, 3, 4], pa.int64()),
                "string": pa.array(["cat", "dog", "egg"], pa.large_string()),
            },
        ],
        [
            "int32-py-both-none-slice",
            ["int32"],
            None,
            [slice(None, None)],
            "default01234",
        ],
        [
            "int32-numpy-untyped",
            ["int32"],
            None,
            [np.asarray([3202, 3203])],
            "default23",
        ],
        [
            "int32-numpy-typed",
            ["int32"],
            None,
            [np.asarray([3202, 3203], np.int32)],
            "default23",
        ],
        [
            "int32-pa-array-typed",
            ["int32"],
            None,
            [pa.array([3202, 3203], pa.int32())],
            "default23",
        ],
        # Index by int16
        [
            "INT16-ALL",
            ["int16"],
            None,
            [],
            "default01234",
        ],
        [
            "int16-all-shaped",
            ["int16"],
            [[1600, 1700]],
            [],
            "default01234",
        ],
        [
            "int16-py-list",
            ["int16"],
            None,
            [[1602, 1603]],
            "default23",
        ],
        [
            "int16-py-tuple",
            ["int16"],
            None,
            [[1602, 1603]],
            "default23",
        ],
        [
            "int16-py-slice",
            ["int16"],
            None,
            [slice(1602, 1603)],
            "default23",
        ],
        [
            "int16-py-left-none-slice",
            ["int16"],
            None,
            [slice(None, 1602)],
            {
                "soma_joinid": pa.array([0, 1, 2], pa.int64()),
                "string": pa.array(["apple", "ball", "cat"], pa.large_string()),
            },
        ],
        [
            "int16-py-right-none-slice",
            ["int16"],
            None,
            [slice(1602, None)],
            {
                "soma_joinid": pa.array([2, 3, 4], pa.int64()),
                "string": pa.array(["cat", "dog", "egg"], pa.large_string()),
            },
        ],
        [
            "int16-py-both-none-slice",
            ["int16"],
            None,
            [slice(None, None)],
            "default01234",
        ],
        [
            "int16-numpy-untyped",
            ["int16"],
            None,
            [np.asarray([1602, 1603])],
            "default23",
        ],
        [
            "int16-numpy-typed",
            ["int16"],
            None,
            [np.asarray([1602, 1603], np.int16)],
            "default23",
        ],
        [
            "int16-pa-array-typed",
            ["int16"],
            None,
            [pa.array([1602, 1603], pa.int16())],
            "default23",
        ],
        # Index by int8
        [
            "INT8-ALL",
            ["int8"],
            None,
            [],
            "default01234",
        ],
        [
            "int8-all-shaped",
            ["int8"],
            [[80, 90]],
            [],
            "default01234",
        ],
        [
            "int8-py-list",
            ["int8"],
            None,
            [[82, 83]],
            "default23",
        ],
        [
            "int8-py-tuple",
            ["int8"],
            None,
            [[82, 83]],
            "default23",
        ],
        [
            "int8-py-slice",
            ["int8"],
            None,
            [slice(82, 83)],
            "default23",
        ],
        [
            "int8-py-left-none-slice",
            ["int8"],
            None,
            [slice(None, 82)],
            {
                "soma_joinid": pa.array([0, 1, 2], pa.int64()),
                "string": pa.array(["apple", "ball", "cat"], pa.large_string()),
            },
        ],
        [
            "int8-py-right-none-slice",
            ["int8"],
            None,
            [slice(82, None)],
            {
                "soma_joinid": pa.array([2, 3, 4], pa.int64()),
                "string": pa.array(["cat", "dog", "egg"], pa.large_string()),
            },
        ],
        [
            "int8-py-both-none-slice",
            ["int8"],
            None,
            [slice(None, None)],
            "default01234",
        ],
        [
            "int8-numpy-untyped",
            ["int8"],
            None,
            [np.asarray([82, 83])],
            "default23",
        ],
        [
            "int8-numpy-typed",
            ["int8"],
            None,
            [np.asarray([82, 83], np.int8)],
            "default23",
        ],
        [
            "int8-pa-array-typed",
            ["int8"],
            None,
            [pa.array([82, 83], pa.int8())],
            "default23",
        ],
        # Index by uint64
        [
            "UINT64-ALL",
            ["uint64"],
            None,
            [],
            "default01234",
        ],
        [
            "uint64-all-shaped",
            ["uint64"],
            [[6410, 6490]],
            [],
            "default01234",
        ],
        [
            "uint64-py-list",
            ["uint64"],
            None,
            [[6412, 6413]],
            "default23",
        ],
        [
            "uint64-py-tuple",
            ["uint64"],
            None,
            [[6412, 6413]],
            "default23",
        ],
        [
            "uint64-py-slice",
            ["uint64"],
            None,
            [slice(6412, 6413)],
            "default23",
        ],
        [
            "uint64-py-left-none-slice",
            ["uint64"],
            None,
            [slice(None, 6412)],
            {
                "soma_joinid": pa.array([0, 1, 2], pa.int64()),
                "string": pa.array(["apple", "ball", "cat"], pa.large_string()),
            },
        ],
        [
            "uint64-py-right-none-slice",
            ["uint64"],
            None,
            [slice(6412, None)],
            {
                "soma_joinid": pa.array([2, 3, 4], pa.int64()),
                "string": pa.array(["cat", "dog", "egg"], pa.large_string()),
            },
        ],
        [
            "uint64-py-both-none-slice",
            ["uint64"],
            None,
            [slice(None, None)],
            "default01234",
        ],
        [
            "uint64-numpy-untyped",
            ["uint64"],
            None,
            [np.asarray([6412, 6413])],
            "default23",
        ],
        [
            "uint64-numpy-typed",
            ["uint64"],
            None,
            [np.asarray([6412, 6413], np.uint64)],
            "default23",
        ],
        [
            "uint64-pa-array-typed",
            ["uint64"],
            None,
            [pa.array([6412, 6413], pa.uint64())],
            "default23",
        ],
        # Index by uint32
        [
            "UINT32-ALL",
            ["uint32"],
            None,
            [],
            "default01234",
        ],
        [
            "uint32-all-shaped",
            ["uint32"],
            [[3210, 3310]],
            [],
            "default01234",
        ],
        [
            "uint32-py-list",
            ["uint32"],
            None,
            [[3212, 3213]],
            "default23",
        ],
        [
            "uint32-py-tuple",
            ["uint32"],
            None,
            [[3212, 3213]],
            "default23",
        ],
        [
            "uint32-py-slice",
            ["uint32"],
            None,
            [slice(3212, 3213)],
            "default23",
        ],
        [
            "uint32-py-left-none-slice",
            ["uint32"],
            None,
            [slice(None, 3212)],
            {
                "soma_joinid": pa.array([0, 1, 2], pa.int64()),
                "string": pa.array(["apple", "ball", "cat"], pa.large_string()),
            },
        ],
        [
            "uint32-py-right-none-slice",
            ["uint32"],
            None,
            [slice(3212, None)],
            {
                "soma_joinid": pa.array([2, 3, 4], pa.int64()),
                "string": pa.array(["cat", "dog", "egg"], pa.large_string()),
            },
        ],
        [
            "uint32-py-both-none-slice",
            ["uint32"],
            None,
            [slice(None, None)],
            "default01234",
        ],
        [
            "uint32-numpy-untyped",
            ["uint32"],
            None,
            [np.asarray([3212, 3213])],
            "default23",
        ],
        [
            "uint32-numpy-typed",
            ["uint32"],
            None,
            [np.asarray([3212, 3213], np.uint32)],
            "default23",
        ],
        [
            "uint32-pa-array-typed",
            ["uint32"],
            None,
            [pa.array([3212, 3213], pa.uint32())],
            "default23",
        ],
        # Index by uint16
        [
            "UINT16-ALL",
            ["uint16"],
            None,
            [],
            "default01234",
        ],
        [
            "uint16-all-shaped",
            ["uint16"],
            [[1610, 1620]],
            [],
            "default01234",
        ],
        [
            "uint16-py-list",
            ["uint16"],
            None,
            [[1612, 1613]],
            "default23",
        ],
        [
            "uint16-py-tuple",
            ["uint16"],
            None,
            [[1612, 1613]],
            "default23",
        ],
        [
            "uint16-py-slice",
            ["uint16"],
            None,
            [slice(1612, 1613)],
            "default23",
        ],
        [
            "uint16-py-left-none-slice",
            ["uint16"],
            None,
            [slice(None, 1612)],
            {
                "soma_joinid": pa.array([0, 1, 2], pa.int64()),
                "string": pa.array(["apple", "ball", "cat"], pa.large_string()),
            },
        ],
        [
            "uint16-py-right-none-slice",
            ["uint16"],
            None,
            [slice(1612, None)],
            {
                "soma_joinid": pa.array([2, 3, 4], pa.int64()),
                "string": pa.array(["cat", "dog", "egg"], pa.large_string()),
            },
        ],
        [
            "uint16-py-both-none-slice",
            ["uint16"],
            None,
            [slice(None, None)],
            "default01234",
        ],
        [
            "uint16-numpy-untyped",
            ["uint16"],
            None,
            [np.asarray([1612, 1613])],
            "default23",
        ],
        [
            "uint16-numpy-typed",
            ["uint16"],
            None,
            [np.asarray([1612, 1613], np.uint16)],
            "default23",
        ],
        [
            "uint16-pa-array-typed",
            ["uint16"],
            None,
            [pa.array([1612, 1613], pa.uint16())],
            "default23",
        ],
        # Index by uint8
        [
            "UINT8-ALL",
            ["uint8"],
            None,
            [],
            "default01234",
        ],
        [
            "uint8-all-shaped",
            ["uint8"],
            [[90, 100]],
            [],
            "default01234",
        ],
        [
            "uint8-py-list",
            ["uint8"],
            None,
            [[92, 93]],
            "default23",
        ],
        [
            "uint8-py-tuple",
            ["uint8"],
            None,
            [[92, 93]],
            "default23",
        ],
        [
            "uint8-py-slice",
            ["uint8"],
            None,
            [slice(92, 93)],
            "default23",
        ],
        [
            "uint8-py-left-none-slice",
            ["uint8"],
            None,
            [slice(None, 92)],
            {
                "soma_joinid": pa.array([0, 1, 2], pa.int64()),
                "string": pa.array(["apple", "ball", "cat"], pa.large_string()),
            },
        ],
        [
            "uint8-py-right-none-slice",
            ["uint8"],
            None,
            [slice(92, None)],
            {
                "soma_joinid": pa.array([2, 3, 4], pa.int64()),
                "string": pa.array(["cat", "dog", "egg"], pa.large_string()),
            },
        ],
        [
            "uint8-py-both-none-slice",
            ["uint8"],
            None,
            [slice(None, None)],
            "default01234",
        ],
        [
            "uint8-numpy-untyped",
            ["uint8"],
            None,
            [np.asarray([92, 93])],
            "default23",
        ],
        [
            "uint8-numpy-typed",
            ["uint8"],
            None,
            [np.asarray([92, 93], np.uint8)],
            "default23",
        ],
        [
            "uint8-pa-array-typed",
            ["uint8"],
            None,
            [pa.array([92, 93], pa.uint8())],
            "default23",
        ],
        # Index by float32
        [
            "FLOAT32-ALL",
            ["float32"],
            None,
            [],
            "default01234",
        ],
        [
            "float32-all-shaped",
            ["float32"],
            [[320.0, 330.0]],
            [],
            "default01234",
        ],
        [
            "float32-py-list",
            ["float32"],
            None,
            [[322.5, 323.5]],
            "default23",
        ],
        [
            "float32-py-tuple",
            ["float32"],
            None,
            [(322.5, 323.5)],
            "default23",
        ],
        [
            "float32-py-slice",
            ["float32"],
            None,
            [slice(322.5, 323.5)],
            "default23",
        ],
        [
            "float32-np-array-untyped",
            ["float32"],
            None,
            [np.asarray([322.5, 323.5])],
            "default23",
        ],
        [
            "float32-np-array-typed",
            ["float32"],
            None,
            [np.asarray([322.5, 323.5], np.float32)],
            "default23",
        ],
        [
            "float32-pa-array-typed-float32",
            ["float32"],
            None,
            [pa.array([322.5, 323.5], pa.float32())],
            "default23",
        ],
        # Index by float64
        [
            "FLOAT64-ALL",
            ["float64"],
            None,
            [],
            "default01234",
        ],
        [
            "float64-all-shaped",
            ["float64"],
            [[640.0, 650.0]],
            [],
            "default01234",
        ],
        [
            "float64-py-list",
            ["float64"],
            None,
            [[642.5, 643.5]],
            "default23",
        ],
        [
            "float64-py-tuple",
            ["float64"],
            None,
            [(642.5, 643.5)],
            "default23",
        ],
        [
            "float64-py-slice",
            ["float64"],
            None,
            [slice(642.5, 643.5)],
            "default23",
        ],
        [
            "float64-np-array-untyped",
            ["float64"],
            None,
            [np.asarray([642.5, 643.5])],
            "default23",
        ],
        [
            "float64-np-array-typed",
            ["float64"],
            None,
            [np.asarray([642.5, 643.5], np.float64)],
            "default23",
        ],
        [
            "float64-pa-array-untyped",
            ["float64"],
            None,
            [pa.array([642.5, 643.5])],
            "default23",
        ],
        [
            "float64-pa-array-typed-float64",
            ["float64"],
            None,
            [pa.array([642.5, 643.5], pa.float64())],
            "default23",
        ],
        # Index by int64 and string
        [
            "INT64+STRING-ALL",
            ["int64", "string"],
            None,
            [],
            "default01234",
        ],
        [
            "int64+string-all-shaped",
            ["int64", "string"],
            [[6400, 6500], None],
            [],
            "default01234",
        ],
        [
            "int64+string-arrow",
            ["int64", "string"],
            None,
            [pa.array([6402, 6403]), pa.array(["cat", "dog"])],
            "default23",
        ],
        [
            "string+int64-arrow",
            ["string", "int64"],
            None,
            [pa.array(["cat", "dog"]), pa.array([6402, 6403])],
            "default23",
        ],
        [
            "string+int64-numpy",
            ["string", "int64"],
            None,
            [np.asarray(["cat", "dog"]), np.asarray([6402, 6403])],
            "default23",
        ],
        [
            "string+int64-py-list",
            ["string", "int64"],
            None,
            [["cat", "dog"], [6402, 6403]],
            "default23",
        ],
        [
            "string+int64-py-tuple",
            ["string", "int64"],
            None,
            [("cat", "dog"), (6402, 6403)],
            "default23",
        ],
        # Index by int64, float64, and string
        [
            "INT64+FLOAT64+STRING-ALL",
            ["int64", "float64", "string"],
            None,
            [],
            "default01234",
        ],
        [
            "int64+float64+string-arrow",
            ["int64", "float64", "string"],
            None,
            [
                pa.array([6402, 6403]),
                pa.array([642.5, 643.5]),
                pa.array(["cat", "dog"]),
            ],
            "default23",
        ],
        [
            "float64+string+int64-arrow",
            ["float64", "string", "int64"],
            None,
            [
                pa.array([642.5, 643.5]),
                pa.array(["cat", "dog"]),
                pa.array([6402, 6403]),
            ],
            "default23",
        ],
        [
            "string+int64+float64-numpy",
            ["string", "int64", "float64"],
            None,
            [
                np.asarray(["cat", "dog"]),
                np.asarray([6402, 6403]),
                np.asarray([642.5, 643.5]),
            ],
            "default23",
        ],
        [
            "string+int64+float64py-list",
            ["string", "int64", "float64"],
            None,
            [["cat", "dog"], [6402, 6403], [642.5, 643.5]],
            "default23",
        ],
        [
            "string+int64+float64-py-tuple",
            ["string", "int64", "float64"],
            None,
            [("cat", "dog"), (6402, 6403), (642.5, 643.5)],
            "default23",
        ],
        # Index by timestamp-s
        [
            "TIMESTAMP-SEC-ALL",
            ["tss"],
            None,
            [],
            "default01234",
        ],
        [
            "timestamp-sec-all-shaped",
            ["tss"],
            [
                [
                    np.datetime64(946684800, "s"),
                    np.datetime64(946684809, "s"),
                ]
            ],
            [],
            "default01234",
        ],
        [
            "tss-py-list",
            ["tss"],
            None,
            [[np.datetime64(946684802, "s"), np.datetime64(946684803, "s")]],
            "default23",
        ],
        [
            "tss-py-tuple",
            ["tss"],
            None,
            [[np.datetime64(946684802, "s"), np.datetime64(946684803, "s")]],
            "default23",
        ],
        [
            "tss-py-slice",
            ["tss"],
            None,
            [slice(np.datetime64(946684802, "s"), np.datetime64(946684803, "s"))],
            "default23",
        ],
        [
            "tss-py-left-none-slice",
            ["tss"],
            None,
            [slice(None, np.datetime64(946684802, "s"))],
            {
                "soma_joinid": pa.array([0, 1, 2], pa.int64()),
                "string": pa.array(["apple", "ball", "cat"], pa.large_string()),
            },
        ],
        [
            "tss-py-right-none-slice",
            ["tss"],
            None,
            [slice(np.datetime64(946684802, "s"), None)],
            {
                "soma_joinid": pa.array([2, 3, 4], pa.int64()),
                "string": pa.array(["cat", "dog", "egg"], pa.large_string()),
            },
        ],
        [
            "tss-py-both-none-slice",
            ["tss"],
            None,
            [slice(None, None)],
            "default01234",
        ],
        [
            "tss-numpy",
            ["tss"],
            None,
            [
                np.asarray(
                    [np.datetime64(946684802, "s"), np.datetime64(946684803, "s")]
                )
            ],
            "default23",
        ],
        [
            "tss-pa-array-untyped",
            ["tss"],
            None,
            [pa.array([946684802, 946684803])],
            "default23",
        ],
        [
            "tss-pa-array-typed",
            ["tss"],
            None,
            [pa.array([946684802, 946684803], pa.timestamp("s"))],
            "default23",
        ],
        # Index by timestamp-ms
        [
            "TIMESTAMP-MSEC-ALL",
            ["tsms"],
            None,
            [],
            "default01234",
        ],
        [
            "timestamp-msec-all-shaped",
            ["tsms"],
            [
                [
                    np.datetime64(946684800000, "ms"),
                    np.datetime64(946684800009, "ms"),
                ]
            ],
            [],
            "default01234",
        ],
        [
            "tsms-py-list",
            ["tsms"],
            None,
            [[np.datetime64(946684800002, "ms"), np.datetime64(946684800003, "ms")]],
            "default23",
        ],
        [
            "tsms-py-tuple",
            ["tsms"],
            None,
            [[np.datetime64(946684800002, "ms"), np.datetime64(946684800003, "ms")]],
            "default23",
        ],
        [
            "tsms-py-slice",
            ["tsms"],
            None,
            [
                slice(
                    np.datetime64(946684800002, "ms"), np.datetime64(946684800003, "ms")
                )
            ],
            "default23",
        ],
        [
            "tsms-py-left-none-slice",
            ["tsms"],
            None,
            [slice(None, np.datetime64(946684800002, "ms"))],
            {
                "soma_joinid": pa.array([0, 1, 2], pa.int64()),
                "string": pa.array(["apple", "ball", "cat"], pa.large_string()),
            },
        ],
        [
            "tsms-py-right-none-slice",
            ["tsms"],
            None,
            [slice(np.datetime64(946684800002, "ms"), None)],
            {
                "soma_joinid": pa.array([2, 3, 4], pa.int64()),
                "string": pa.array(["cat", "dog", "egg"], pa.large_string()),
            },
        ],
        [
            "tsms-py-both-none-slice",
            ["tsms"],
            None,
            [slice(None, None)],
            "default01234",
        ],
        [
            "tsms-numpy",
            ["tsms"],
            None,
            [
                np.asarray(
                    [
                        np.datetime64(946684800002, "ms"),
                        np.datetime64(946684800003, "ms"),
                    ]
                )
            ],
            "default23",
        ],
        [
            "tsms-pa-array-untyped",
            ["tsms"],
            None,
            [pa.array([946684800002, 946684800003])],
            "default23",
        ],
        [
            "tsms-pa-array-typed",
            ["tsms"],
            None,
            [pa.array([946684800002, 946684800003], pa.timestamp("ms"))],
            "default23",
        ],
        # Index by timestamp-us
        [
            "TIMESTAMP-USEC-ALL",
            ["tsus"],
            None,
            [],
            "default01234",
        ],
        [
            "timestamp-usec-all-shaped",
            ["tsus"],
            [
                [
                    np.datetime64(946684800000000, "us"),
                    np.datetime64(946684800000009, "us"),
                ]
            ],
            [],
            "default01234",
        ],
        [
            "tsus-py-list",
            ["tsus"],
            None,
            [
                [
                    np.datetime64(946684800000002, "us"),
                    np.datetime64(946684800000003, "us"),
                ]
            ],
            "default23",
        ],
        [
            "tsus-py-tuple",
            ["tsus"],
            None,
            [
                [
                    np.datetime64(946684800000002, "us"),
                    np.datetime64(946684800000003, "us"),
                ]
            ],
            "default23",
        ],
        [
            "tsus-py-slice",
            ["tsus"],
            None,
            [
                slice(
                    np.datetime64(946684800000002, "us"),
                    np.datetime64(946684800000003, "us"),
                )
            ],
            "default23",
        ],
        [
            "tsus-py-left-none-slice",
            ["tsus"],
            None,
            [slice(None, np.datetime64(946684800000002, "us"))],
            {
                "soma_joinid": pa.array([0, 1, 2], pa.int64()),
                "string": pa.array(["apple", "ball", "cat"], pa.large_string()),
            },
        ],
        [
            "tsus-py-right-none-slice",
            ["tsus"],
            None,
            [slice(np.datetime64(946684800000002, "us"), None)],
            {
                "soma_joinid": pa.array([2, 3, 4], pa.int64()),
                "string": pa.array(["cat", "dog", "egg"], pa.large_string()),
            },
        ],
        [
            "tsus-py-both-none-slice",
            ["tsus"],
            None,
            [slice(None, None)],
            "default01234",
        ],
        [
            "tsus-numpy",
            ["tsus"],
            None,
            [
                np.asarray(
                    [
                        np.datetime64(946684800000002, "us"),
                        np.datetime64(946684800000003, "us"),
                    ]
                )
            ],
            "default23",
        ],
        [
            "tsus-pa-array-untyped",
            ["tsus"],
            None,
            [pa.array([946684800000002, 946684800000003])],
            "default23",
        ],
        [
            "tsus-pa-array-typed",
            ["tsus"],
            None,
            [pa.array([946684800000002, 946684800000003], pa.timestamp("us"))],
            "default23",
        ],
        # Index by timestamp-ns
        [
            "TIMESTAMP-NSEC-ALL",
            ["tsns"],
            None,
            [],
            "default01234",
        ],
        [
            "timestamp-nsec-all-shaped",
            ["tsns"],
            [
                [
                    np.datetime64(946684800000000000, "ns"),
                    np.datetime64(946684800000000009, "ns"),
                ]
            ],
            [],
            "default01234",
        ],
        [
            "tsns-py-list",
            ["tsns"],
            None,
            [
                [
                    np.datetime64(946684800000000002, "ns"),
                    np.datetime64(946684800000000003, "ns"),
                ]
            ],
            "default23",
        ],
        [
            "tsns-py-tuple",
            ["tsns"],
            None,
            [
                [
                    np.datetime64(946684800000000002, "ns"),
                    np.datetime64(946684800000000003, "ns"),
                ]
            ],
            "default23",
        ],
        [
            "tsns-py-slice",
            ["tsns"],
            None,
            [
                slice(
                    np.datetime64(946684800000000002, "ns"),
                    np.datetime64(946684800000000003, "ns"),
                )
            ],
            "default23",
        ],
        [
            "tsns-py-left-none-slice",
            ["tsns"],
            None,
            [slice(None, np.datetime64(946684800000000002, "ns"))],
            {
                "soma_joinid": pa.array([0, 1, 2], pa.int64()),
                "string": pa.array(["apple", "ball", "cat"], pa.large_string()),
            },
        ],
        [
            "tsns-py-right-none-slice",
            ["tsns"],
            None,
            [slice(np.datetime64(946684800000000002, "ns"), None)],
            {
                "soma_joinid": pa.array([2, 3, 4], pa.int64()),
                "string": pa.array(["cat", "dog", "egg"], pa.large_string()),
            },
        ],
        [
            "tsns-py-both-none-slice",
            ["tsns"],
            None,
            [slice(None, None)],
            "default01234",
        ],
        [
            "tsns-numpy",
            ["tsns"],
            None,
            [
                np.asarray(
                    [
                        np.datetime64(946684800000000002, "ns"),
                        np.datetime64(946684800000000003, "ns"),
                    ]
                )
            ],
            "default23",
        ],
        [
            "tsns-pa-array-untyped",
            ["tsns"],
            None,
            [pa.array([946684800000000002, 946684800000000003])],
            "default23",
        ],
        [
            "tsns-pa-array-typed",
            ["tsns"],
            None,
            [pa.array([946684800000000002, 946684800000000003], pa.timestamp("ns"))],
            "default23",
        ],
    ],
)
def test_types_no_errors(
    tmp_path,
    arrow_table,
    name,
    index_column_names,
    domain,
    coords,
    expecteds,
):
    uri = tmp_path.as_posix()

    soma.DataFrame.create(
        uri,
        schema=arrow_table.schema,
        index_column_names=index_column_names,
        domain=domain,
    )

    with soma.DataFrame.open(uri, "w") as sdf:
        sdf.write(arrow_table)

    if expecteds == "default01234":
        expecteds = {
            "soma_joinid": pa.array([0, 1, 2, 3, 4], pa.int64()),
            "string": pa.array(
                ["apple", "ball", "cat", "dog", "egg"], pa.large_string()
            ),
        }

    elif expecteds == "default23":
        expecteds = {
            "soma_joinid": pa.array([2, 3], pa.int64()),
            "string": pa.array(["cat", "dog"], pa.large_string()),
        }

    with soma.DataFrame.open(uri, "r") as sdf:
        actual_table = sdf.read(coords=coords).concat()
        for query_column_name, expected_array in expecteds.items():
            actual_array = actual_table[query_column_name].combine_chunks()
            # The output from the first one is easier to read when it fails
            assert actual_array.to_pylist() == expected_array.to_pylist()
            assert actual_array == expected_array
        if domain is not None:
            actual_domain = sdf.domain
            for i in range(len(domain)):
                if domain[i] is not None:
                    assert actual_domain[i] == tuple(domain[i])
        else:
            # TileDB always returns ("", "") as domain for string/bytes dimensions on sparse arrays.
            actual_domain = sdf.domain
            for i in range(len(index_column_names)):
                if (
                    index_column_names[i] == "string"
                    or index_column_names[i] == "bytes"
                ):
                    assert actual_domain[i] == ("", "")


@pytest.mark.parametrize(
    "name,index_column_names,domain,coords,error",
    [
        [
            "soma-joinid-py-list-shaped-negatively",
            ["soma_joinid"],
            [[-100, 100]],
            [],
            ValueError,
        ],
        [
            "soma_joinid-all-shaped-too-short",
            ["soma_joinid"],
            [],
            [],
            ValueError,
        ],
        [
            "soma_joinid-all-shaped-too-long",
            ["soma_joinid"],
            [[0, 10], [0, 10]],
            [],
            ValueError,
        ],
        [
            "string-py-list-shaped-at-all",
            ["string"],
            [["a", "z"]],
            [],
            ValueError,
        ],
        [
            "BOOL-ALL",
            ["bool"],
            None,
            [],
            TypeError,
        ],
    ],
)
def test_types_create_errors(
    tmp_path,
    arrow_table,
    name,
    index_column_names,
    domain,
    coords,
    error,
):
    uri = tmp_path.as_posix()

    with pytest.raises(error):
        soma.DataFrame.create(
            uri,
            schema=arrow_table.schema,
            index_column_names=index_column_names,
            domain=domain,
        )


@pytest.mark.parametrize(
    "name,index_column_names,domain,error",
    [
        [
            "int32-py-list-shaped-out-of-bounds",
            ["int32"],
            [[100, 200]],
            soma.SOMAError,
        ],
        [
            "int16-py-list-shaped-out-of-bounds",
            ["int16"],
            [[100, 200]],
            soma.SOMAError,
        ],
        [
            "int8-py-list-shaped-out-of-bounds",
            ["int8"],
            [[10, 20]],
            soma.SOMAError,
        ],
        [
            "uint64-py-list-shaped-out-of-bounds",
            ["uint64"],
            [[100, 200]],
            soma.SOMAError,
        ],
        [
            "uint32-py-list-shaped-out-of-bounds",
            ["uint32"],
            [[100, 200]],
            soma.SOMAError,
        ],
        [
            "uint32-py-list-shaped-out-of-bounds",
            ["uint32"],
            [[100, 200]],
            soma.SOMAError,
        ],
        [
            "uint8-py-list-shaped-out-of-bounds",
            ["uint8"],
            [[10, 20]],
            soma.SOMAError,
        ],
        [
            "float32-py-list-shaped-out-of-bounds",
            ["float32"],
            [[100.0, 200.0]],
            soma.SOMAError,
        ],
        [
            "float64-py-list-shaped-out-of-bounds",
            ["float64"],
            [[100.0, 200.0]],
            soma.SOMAError,
        ],
    ],
)
def test_types_write_errors(
    tmp_path,
    arrow_table,
    name,
    index_column_names,
    domain,
    error,
):
    uri = tmp_path.as_posix()

    soma.DataFrame.create(
        uri,
        schema=arrow_table.schema,
        index_column_names=index_column_names,
        domain=domain,
    )

    with pytest.raises(error):
        with soma.DataFrame.open(uri, "w") as sdf:
            sdf.write(arrow_table)


@pytest.mark.parametrize(
    "name,index_column_names,domain,coords",
    [
        [
            "int32-pa-array-untyped",
            ["int32"],
            None,
            [pa.array([3202, 3203])],
            # Static type (INT64) does not match expected type (INT32)
        ],
        [
            "int16-pa-array-untyped",
            ["int16"],
            None,
            [pa.array([1602, 1603])],
            # Static type (INT64) does not match expected type (INT16)
        ],
        [
            "int8-pa-array-untyped",
            ["int8"],
            None,
            [pa.array([82, 83])],
            # Static type (INT64) does not match expected type (INT8)
        ],
        [
            "uint64-pa-array-untyped",
            ["uint64"],
            None,
            [pa.array([6412, 6413])],
            # Static type (INT64) does not match expected type (UINT64)
        ],
        [
            "uint32-pa-array-untyped",
            ["uint32"],
            None,
            [pa.array([3212, 3213])],
            # Static type (UINT64) does not match expected type (UINT32)
        ],
        [
            "uint16-pa-array-untyped",
            ["uint16"],
            None,
            [pa.array([1612, 1613])],
            # Static type (UINT64) does not match expected type (UINT16)
        ],
        [
            "uint8-pa-array-untyped",
            ["uint8"],
            None,
            [pa.array([92, 93])],
            # Static type (UINT64) does not match expected type (UINT8)
        ],
        [
            "float32-pa-array-untyped",
            ["float32"],
            None,
            [pa.array([322.5, 323.5])],
            # Static type (FLOAT64) does not match expected type (FLOAT32)
        ],
        [
            "float32-pa-array-typed-float64",
            ["float32"],
            None,
            [pa.array([322.5, 323.5], pa.float64())],
            # Static type (FLOAT64) does not match expected type (FLOAT32)
        ],
        [
            "float64-pa-array-typed-float32",
            ["float64"],
            None,
            [pa.array([322.5, 323.5], pa.float32())],
            # Static type (FLOAT32) does not match expected type (FLOAT64)
        ],
    ],
)
def test_types_read_errors(
    tmp_path,
    arrow_table,
    name,
    index_column_names,
    domain,
    coords,
):
    uri = tmp_path.as_posix()

    soma.DataFrame.create(
        uri,
        schema=arrow_table.schema,
        index_column_names=index_column_names,
        domain=domain,
    )

    with soma.DataFrame.open(uri, "w") as sdf:
        sdf.write(arrow_table)

    with pytest.raises((soma.SOMAError)):
        with soma.DataFrame.open(uri, "r") as sdf:
            sdf.read(coords=coords).concat()
