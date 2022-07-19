import pyarrow as pa
import pytiledbsc
import pytest
import numpy as np

DATA_SIZE = 8
VERBOSE = False

rng = np.random.default_rng()


def test_int32():
    data = np.random.randint(-1 << 31, 1 << 31, size=DATA_SIZE, dtype=np.int32)
    cb = pytiledbsc.ColumnBuffer("buf", pytiledbsc.DataType.INT32, len(data), data)

    if VERBOSE:
        print(f"Checking: {data} == {cb.data()}")
    assert np.array_equal(data, cb.data())


def test_int64():
    data = np.random.randint(-1 << 63, 1 << 63, size=DATA_SIZE, dtype=np.int64)
    cb = pytiledbsc.ColumnBuffer("buf", pytiledbsc.DataType.INT64, len(data), data)

    if VERBOSE:
        print(f"Checking: {data} == {cb.data()}")
    assert np.array_equal(data, cb.data())


def test_float32():
    data = rng.random(size=DATA_SIZE, dtype=np.float32)
    cb = pytiledbsc.ColumnBuffer("buf", pytiledbsc.DataType.FLOAT32, len(data), data)

    if VERBOSE:
        print(f"Checking: {data} == {cb.data()}")
    assert np.array_equal(data, cb.data())


def test_float64():
    data = rng.random(size=DATA_SIZE, dtype=np.float64)
    cb = pytiledbsc.ColumnBuffer("buf", pytiledbsc.DataType.FLOAT64, len(data), data)

    if VERBOSE:
        print(f"Checking: {data} == {cb.data()}")
    assert np.array_equal(data, cb.data())


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

        if VERBOSE:
            print(f"Checking: {buf.data()} == {data}")
            print(f"Checking: {buf.offsets()} == {offsets}")
            print(f"Checking: {buf.validity()} == {validity}")

        assert np.array_equal(buf.data(), data)
        assert np.array_equal(buf.offsets(), offsets)
        assert np.array_equal(buf.validity(), validity)


def cb_to_arrow(cb):
    return pa.Array._import_from_c(*cb.to_arrow())


def test_arrow():
    data = np.random.randint(-1 << 31, 1 << 31, size=DATA_SIZE, dtype=np.int32)
    cb = pytiledbsc.ColumnBuffer("buf", pytiledbsc.DataType.INT32, len(data), data)

    if VERBOSE:
        print(f"Checking: {data} == {cb.data()}")
    assert np.array_equal(data, cb.data())

    arrow = cb_to_arrow(cb)
    assert np.array_equal(data, arrow)


if __name__ == "__main__":
    test_arrow()
