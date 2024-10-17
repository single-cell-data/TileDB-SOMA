import pathlib

import pyarrow as pa
import pytest

import tiledbsoma as soma

"""
Arrow types we expect to work. A handful of types will promote, e.g., string->large_string.
Most must be literally as requested, or error out.

Tuple is (requested_type, expected_type).
"""
SUPPORTED_ARROW_TYPES = [
    (pa.bool_(),) * 2,
    (pa.int8(),) * 2,
    (pa.uint8(),) * 2,
    (pa.int16(),) * 2,
    (pa.uint16(),) * 2,
    (pa.int32(),) * 2,
    (pa.uint32(),) * 2,
    (pa.int64(),) * 2,
    (pa.uint64(),) * 2,
    (pa.float32(),) * 2,
    (pa.float64(),) * 2,
    (pa.timestamp("s"),) * 2,
    (pa.timestamp("ms"),) * 2,
    (pa.timestamp("us"),) * 2,
    (pa.timestamp("ns"),) * 2,
    (pa.string(), pa.large_string()),
    (pa.binary(), pa.large_binary()),
    (pa.large_string(),) * 2,
    (pa.large_binary(),) * 2,
    (pa.dictionary(pa.int32(), pa.string()),) * 2,
]


"""Arrow types we expect to fail"""
UNSUPPORTED_ARROW_TYPES = [
    pa.null(),
    pa.date64(),
    pa.time64("us"),
    pa.time64("ns"),
    pa.float16(),
    pa.date32(),
    pa.time32("s"),
    pa.time32("ms"),
    pa.duration("s"),
    pa.duration("ms"),
    pa.duration("us"),
    pa.duration("ns"),
    pa.month_day_nano_interval(),
    pa.binary(10),
    pa.decimal128(1),
    pa.decimal128(38),
    pa.list_(pa.int8()),
    pa.large_list(pa.bool_()),
    pa.map_(pa.string(), pa.int32()),
    pa.struct([("f1", pa.int32()), ("f2", pa.string())]),
]


@pytest.mark.parametrize("arrow_type_info", SUPPORTED_ARROW_TYPES)
def test_arrow_types_supported(tmp_path: pathlib.Path, arrow_type_info):
    """Verify round-trip conversion of types which should work "as is" """
    arrow_type, expected_arrow_type = arrow_type_info
    sdf = soma.DataFrame.create(
        tmp_path.as_posix(), schema=pa.schema([(str(arrow_type), arrow_type)])
    )
    schema = sdf.schema
    assert schema is not None
    assert sorted(schema.names) == sorted(["soma_joinid", str(arrow_type)])
    assert schema.field(str(arrow_type)).type == expected_arrow_type


@pytest.mark.parametrize("arrow_type", UNSUPPORTED_ARROW_TYPES)
def test_arrow_types_unsupported(tmp_path, arrow_type):
    """Verify explicit error for unsupported types"""

    with pytest.raises(TypeError, match=r"unsupported type|Unsupported Arrow type"):
        soma.DataFrame.create(
            tmp_path.as_posix(), schema=pa.schema([(str(arrow_type), arrow_type)])
        )


# ================================================================
# Test writing bool byte-arrays to Arrow bit-arrays
@pytest.mark.parametrize(
    "bool_array",
    [
        # Length-zero bit-array
        [],
        # Length-one bit-array
        [True],
        # Less than a full byte
        [True, False, False, True, False, False, False],
        # A single byte
        [True, False, False, True, False, False, False, True],
        # More than a single byte
        [True, False, False, True, False, False, False, True, True],
        # Multiple bytes
        [
            True,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            True,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            True,
            False,
            False,
            False,
            False,
            False,
        ],
    ],
)
def test_bool_arrays(tmp_path, bool_array):
    schema = pa.schema(
        [
            ("soma_joinid", pa.int64()),
            ("b", pa.bool_()),
        ]
    )
    index_column_names = ["soma_joinid"]

    n_data = len(bool_array)
    data = {
        "soma_joinid": list(range(n_data)),
        "b": bool_array,
    }
    rb = pa.Table.from_pydict(data)
    nrb = len(rb)

    with soma.DataFrame.create(
        tmp_path.as_posix(),
        schema=schema,
        index_column_names=index_column_names,
        domain=[[0, max(nrb, 1) - 1]],
    ) as sdf:
        sdf.write(rb)

    with soma.DataFrame.open(tmp_path.as_posix()) as sdf:
        table = sdf.read().concat()
    assert table["b"].to_pylist() == bool_array
