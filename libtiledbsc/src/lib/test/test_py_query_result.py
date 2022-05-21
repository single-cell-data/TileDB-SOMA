import numpy as np
import pyarrow as pa

import pytest
import py_query_result_aux as qra

from typing import Union

from numpy.testing import assert_array_equal


def compare_maybe_ndarray(a: Union[None, np.ndarray], b: Union[None, np.ndarray]):
    if not (isinstance(a, np.ndarray) and isinstance(b, np.ndarray)):
        return a == b
    return np.array_equal(a, b)


def test_bufferset_init():
    orig_data = np.array([1, 2, 3, 4], dtype=np.int32)
    orig_offsets = np.array([0, 0, 2, 3, 3], dtype=np.uint64)
    orig_validity = np.array([1, 0, 0, 1], dtype=np.uint8)
    ov_pairs = [
        (orig_data, None, None),
        (orig_data, orig_offsets, None),
        (orig_data, orig_offsets, orig_validity),
    ]

    for ov in ov_pairs:
        data, offsets, validity = ov
        buf = qra.BufferSet("b", qra.DataType.INT64, 8, data, offsets, validity)

        assert_array_equal(buf.data(), data.view(np.uint8))
        assert compare_maybe_ndarray(buf.offsets(), offsets)
        assert compare_maybe_ndarray(buf.validity(), validity)

    # elem_nbytes does not match size of datatype
    # TODO elem_nbytes should actually be cell_num
    with pytest.raises(Exception):
        # TODO assert TileDBSCError
        qra.BufferSet("b", qra.DataType.INT64, 4, data)


def test_query_result_basic():
    orig_data = np.array([1, 2, 3, 4], dtype=np.int32)
    orig_validity = np.array([1, 0, 0, 1], dtype=np.uint8)

    buf = qra.BufferSet("b", qra.DataType.INT32, 4, orig_data, None, orig_validity)

    res = qra.QueryResult({"b": buf})

    arw = res.to_arrow("b")
    assert isinstance(arw, qra.ArrowPair)
    arrow_array = pa.Array._import_from_c(arw.array(), arw.schema())
    # TODO should this comparison use None instead of np.nan?
    assert_array_equal(arrow_array, [1, np.nan, np.nan, 4])


def test_query_result_to_arrow_table():
    def create_bufferset(name, data, type, offsets, validity):
        return qra.BufferSet(name, type, 1, data, offsets, None)

    validity = np.ones(11, dtype=np.uint8)
    validity_idxs = [1, 3, 5, 6, 7]
    validity[validity_idxs] = 0

    strings = ["", "abcd", "ef", "ghijk", "lmno", "p", "", "q", "rstu", "vwxyz", ""]
    pa_strings = pa.array(strings)
    (
        str_offsets,
        str_data,
    ) = map(np.array, pa_strings.buffers()[1:])
    str_offsets = str_offsets.view(np.uint32).astype(np.uint64)
    strings_nulled = np.array(strings, dtype="O")
    strings_nulled[validity_idxs] = None

    int_data = np.arange(11, dtype=np.int32)
    int_data_nulled = np.array(int_data, dtype=np.float32)
    int_data_nulled[validity_idxs] = np.nan

    res = qra.QueryResult(
        {
            "ints": qra.BufferSet("ints", qra.DataType.INT32, 4, int_data, None, None),
            "nullable_ints": qra.BufferSet(
                "nullable_ints", qra.DataType.INT32, 4, int_data, None, validity
            ),
            "strings": qra.BufferSet(
                "strings", qra.DataType.STRING_ASCII, 1, str_data, str_offsets, None
            ),
            "nullable_strings": qra.BufferSet(
                "nullable_strings",
                qra.DataType.STRING_ASCII,
                1,
                str_data,
                str_offsets,
                validity,
            ),
        }
    )

    arw_ints = res.to_arrow("ints")
    ret_ints = pa.Array._import_from_c(arw_ints.array(), arw_ints.schema())
    assert_array_equal(ret_ints, int_data)

    arw_nl_ints = res.to_arrow("nullable_ints")
    ret_nl_ints = pa.Array._import_from_c(arw_nl_ints.array(), arw_nl_ints.schema())
    assert_array_equal(ret_nl_ints, int_data_nulled)

    arw_strings = res.to_arrow("strings")
    ret_strings = pa.Array._import_from_c(arw_strings.array(), arw_strings.schema())
    assert np.array_equal(strings, ret_strings)

    arw_nl_strings = res.to_arrow("nullable_strings")
    ret_nl_strings = pa.Array._import_from_c(
        arw_nl_strings.array(), arw_nl_strings.schema()
    )
    assert np.array_equal(strings_nulled, ret_nl_strings)
