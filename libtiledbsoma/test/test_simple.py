import random

import numpy as np
import pandas as pd
import pyarrow as pa

import tiledbsoma.pytiledbsoma as clib

DATA_SIZE = 1 << 14
VERBOSE = False

rng = np.random.default_rng()


def random_strings(max_length=32):
    # Generate list of random length strings (omit 0 to avoid string comparison failure)
    chars = "".join(chr(i) for i in range(64, 128))
    strings = [
        "".join(random.choices(chars, k=np.random.randint(max_length)))
        for i in range(DATA_SIZE)
    ]

    # Convert to data and offsets
    pa_data = pa.array(strings)
    offsets, data = map(np.array, pa_data.buffers()[1:])
    offsets = offsets.view(np.uint32).astype(np.uint64)

    cb = clib.ColumnBuffer(
        "buf", clib.DataType.STRING_ASCII, len(strings), data, offsets
    )

    return strings, cb


def check_array(data, cb):
    cb_data = cb.data()

    if VERBOSE:
        print(f"Expected {type(data)} = {data}")
        print(f"ColumnBuffer {type(cb_data)} = {cb_data}")
    assert np.array_equal(data, cb.data())

    arrow = cb.to_arrow()
    if VERBOSE:
        print(f"Arrow: {type(arrow)} = {arrow}")
    assert np.array_equal(data, arrow)


def skip_test_init():
    orig_data = np.random.randint(-1 << 31, 1 << 31, size=DATA_SIZE, dtype=np.int32)
    orig_offsets = np.random.randint(1 << 63, size=DATA_SIZE + 1, dtype=np.uint64)
    orig_validity = np.random.randint(1 << 7, size=DATA_SIZE + 1, dtype=np.uint8)
    ovs = [
        (orig_data, None, None),
        (orig_data, orig_offsets, None),
        (orig_data, orig_offsets, orig_validity),
    ]

    for ov in ovs:
        data, offsets, validity = ov
        buf = clib.ColumnBuffer(
            "buf", clib.DataType.INT32, len(data), data, offsets, validity
        )

        assert np.array_equal(buf.data(), data)
        assert np.array_equal(buf.offsets(), offsets)
        assert np.array_equal(buf.validity(), validity)


def skip_test_int32():
    data = np.random.randint(-1 << 31, 1 << 31, size=DATA_SIZE, dtype=np.int32)
    cb = clib.ColumnBuffer("buf", clib.DataType.INT32, len(data), data)

    check_array(data, cb)


def skip_test_int64():
    data = np.random.randint(-1 << 63, 1 << 63, size=DATA_SIZE, dtype=np.int64)
    cb = clib.ColumnBuffer("buf", clib.DataType.INT64, len(data), data)

    check_array(data, cb)


def skip_test_float32():
    data = rng.random(size=DATA_SIZE, dtype=np.float32)
    cb = clib.ColumnBuffer("buf", clib.DataType.FLOAT32, len(data), data)

    check_array(data, cb)


def skip_test_float64():
    data = rng.random(size=DATA_SIZE, dtype=np.float64)
    cb = clib.ColumnBuffer("buf", clib.DataType.FLOAT64, len(data), data)

    check_array(data, cb)


def skip_test_string():
    strings, cb = random_strings()

    if VERBOSE:
        print(f"Expected {type(strings)} = {strings}")

    arrow = cb.to_arrow()
    if VERBOSE:
        print(f"Arrow: {type(arrow)} = {arrow}")
    assert np.array_equal(strings, arrow)


def skip_test_table():
    a = np.random.randint(-1 << 31, 1 << 31, size=DATA_SIZE, dtype=np.int32)
    a_cb = clib.ColumnBuffer("a", clib.DataType.INT32, len(a), a)
    b = np.random.randint(-1 << 63, 1 << 63, size=DATA_SIZE, dtype=np.int64)
    b_cb = clib.ColumnBuffer("b", clib.DataType.INT64, len(b), b)
    c = rng.random(size=DATA_SIZE, dtype=np.float32)
    c_cb = clib.ColumnBuffer("c", clib.DataType.FLOAT32, len(c), c)
    d = rng.random(size=DATA_SIZE, dtype=np.float64)
    d_cb = clib.ColumnBuffer("d", clib.DataType.FLOAT64, len(d), d)
    e, e_cb = random_strings(8)

    arrow = clib.to_arrow({"a": a_cb, "b": b_cb, "c": c_cb, "d": d_cb, "e": e_cb})
    df = arrow.to_pandas()
    df_expected = pd.DataFrame({"a": a, "b": b, "c": c, "d": d, "e": e})

    if VERBOSE:
        print("Arrow")
        print(type(arrow))
        print(arrow)
        print("Expected")
        print(type(df_expected))
        print(df_expected)
        print("Actual")
        print(type(df))
        print(df)

    assert df_expected.equals(df)


if __name__ == "__main__":
    skip_test_init()
