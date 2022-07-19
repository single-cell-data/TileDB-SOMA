import pyarrow as pa
import pytiledbsc
import pytest
import numpy as np
import random

DATA_SIZE = 8
VERBOSE = True

rng = np.random.default_rng()


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


def test_init():
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
        buf = pytiledbsc.ColumnBuffer(
            "buf", pytiledbsc.DataType.INT32, len(data), data, offsets, validity
        )

        assert np.array_equal(buf.data(), data)
        assert np.array_equal(buf.offsets(), offsets)
        assert np.array_equal(buf.validity(), validity)


def test_int32():
    data = np.random.randint(-1 << 31, 1 << 31, size=DATA_SIZE, dtype=np.int32)
    cb = pytiledbsc.ColumnBuffer("buf", pytiledbsc.DataType.INT32, len(data), data)

    check_array(data, cb)


def test_int64():
    data = np.random.randint(-1 << 63, 1 << 63, size=DATA_SIZE, dtype=np.int64)
    cb = pytiledbsc.ColumnBuffer("buf", pytiledbsc.DataType.INT64, len(data), data)

    check_array(data, cb)


def test_float32():
    data = rng.random(size=DATA_SIZE, dtype=np.float32)
    cb = pytiledbsc.ColumnBuffer("buf", pytiledbsc.DataType.FLOAT32, len(data), data)

    check_array(data, cb)


def test_float64():
    data = rng.random(size=DATA_SIZE, dtype=np.float64)
    cb = pytiledbsc.ColumnBuffer("buf", pytiledbsc.DataType.FLOAT64, len(data), data)

    check_array(data, cb)


def test_string():
    # Generate list of random length strings (omit 0 to avoid string comparison failure)
    max_len = 64
    chars = "".join(chr(i) for i in range(1, 256))
    strings = [
        "".join(random.choices(chars, k=np.random.randint(max_len)))
        for i in range(DATA_SIZE)
    ]

    # Convert to data and offsets
    pa_data = pa.array(strings)
    offsets, data = map(np.array, pa_data.buffers()[1:])
    offsets = offsets.view(np.uint32).astype(np.uint64)

    cb = pytiledbsc.ColumnBuffer(
        "buf", pytiledbsc.DataType.STRING_ASCII, len(strings), data, offsets
    )

    if VERBOSE:
        print(f"Expected {type(strings)} = {strings}")
        print(f"ColumnBuffer {type(data)} = {data}")
        print(f"ColumnBuffer {type(offsets)} = {offsets}")

    arrow = cb.to_arrow()
    if VERBOSE:
        print(f"Arrow: {type(arrow)} = {arrow}")
    assert np.array_equal(strings, arrow)


if __name__ == "__main__":
    test_new()
