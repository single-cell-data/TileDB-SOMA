import numpy as np
import pytiledbsc

DATA_SIZE = 8
rng = np.random.default_rng()


def test_int32():
    data = np.random.randint(-1 << 31, 1 << 31, size=DATA_SIZE, dtype=np.int32)
    cb = pytiledbsc.ColumnBuffer("buf", pytiledbsc.DataType.INT32, len(data), data)

    print(f"Checking: {data} == {cb.data()}")
    assert np.array_equal(data, cb.data())


def test_int64():
    data = np.random.randint(-1 << 63, 1 << 63, size=DATA_SIZE, dtype=np.int64)
    cb = pytiledbsc.ColumnBuffer("buf", pytiledbsc.DataType.INT64, len(data), data)

    print(f"Checking: {data} == {cb.data()}")
    assert np.array_equal(data, cb.data())


def test_float32():
    data = rng.random(size=DATA_SIZE, dtype=np.float32)
    cb = pytiledbsc.ColumnBuffer("buf", pytiledbsc.DataType.FLOAT32, len(data), data)

    print(f"Checking: {data} == {cb.data()}")
    assert np.array_equal(data, cb.data())


def test_float64():
    data = rng.random(size=DATA_SIZE, dtype=np.float64)
    cb = pytiledbsc.ColumnBuffer("buf", pytiledbsc.DataType.FLOAT64, len(data), data)

    print(f"Checking: {data} == {cb.data()}")
    assert np.array_equal(data, cb.data())
