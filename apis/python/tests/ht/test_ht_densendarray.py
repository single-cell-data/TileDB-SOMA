"""Hypothesis tests for DenseNDArray."""

from __future__ import annotations

from typing import Any

import hypothesis.extra.numpy as ht_np
import numpy as np
import pyarrow as pa
import pytest
from hypothesis import strategies as st
from hypothesis.stateful import (
    initialize,
    invariant,
    precondition,
    rule,
)

import tiledbsoma as soma
from tiledbsoma._core_options import OpenMode

from tests.ht._array_state_machine import SOMANDArrayStateMachine
from tests.ht._ht_test_config import HT_TEST_CONFIG
from tests.ht._ht_util import ndarray_datatype
from tests.ht._ledger import ArrowTensorLedgerEntry, Ledger, get_entries


@st.composite
def dense_array_shape(
    draw: st.DrawFn,
    *,
    min_shape: tuple[int, ...] | None = None,
    max_shape: tuple[int, ...] | None = None,
) -> tuple[int | None, ...]:
    """Strategy to generate nd array shapes."""

    MAX_DIMS = 3
    if min_shape is not None:
        ndim = len(min_shape)
    elif max_shape is not None:
        ndim = len(max_shape)
    else:
        ndim = draw(st.integers(min_value=1, max_value=MAX_DIMS))

    min_values = [1] * ndim if min_shape is None else min_shape
    # to keep the array under some number of elements max, use the nth root as max per-dim size
    MAX_ELEM = 2**20  # 1M
    shape_limit = int(MAX_ELEM ** (1 / ndim))
    max_values = [shape_limit] * ndim if max_shape is None else [min(shape_limit, s) for s in max_shape]

    elements = [
        draw(
            (
                st.integers(min_value=min_values[i], max_value=max_values[i])
                if min_values[i] < max_values[i]
                else st.just(min_values[i])
            ),
        )
        for i in range(ndim)
    ]
    new_shape = tuple(elements)

    if min_shape is not None:
        assert len(new_shape) == len(min_shape)
        assert new_shape >= tuple((s or 1) for s in min_shape)

    return new_shape


@st.composite
def dense_indices(draw: st.DrawFn, shape: tuple[int, ...]) -> tuple[int | slice, ...]:
    """Strategy to return DenseNDArray slicing, which currently allows:
    * None - synonym for slice(None)
    * slice - with step == 1 ONLY
    * int - a single integer coord
    """

    def one_dim(s: int) -> int | slice:
        if draw(st.booleans()):
            return draw(st.integers(min_value=0, max_value=s))
        element = st.integers(min_value=0, max_value=s)
        a, b = draw(element), draw(element)
        a, b = (a, b) if a <= b else (b, a)
        if a == 0 and b == s and draw(st.booleans()):
            return slice(None)
        return slice(a, b, None)

    return tuple(one_dim(s) for s in shape)


DEFAULT_FILL_VALUE = {
    pa.int8(): -127,
    pa.int16(): -32767,
    pa.int32(): -2147483647,
    pa.int64(): -9223372036854775807,
    pa.uint8(): 2**8 - 1,
    pa.uint16(): 2**16 - 1,
    pa.uint32(): 2**32 - 1,
    pa.uint64(): 2**64 - 1,
    pa.float32(): np.nan,
    pa.float64(): np.nan,
    pa.bool_(): False,
    pa.timestamp("s"): "NaT",
    pa.timestamp("ms"): "NaT",
    pa.timestamp("us"): "NaT",
    pa.timestamp("ns"): "NaT",
}


def fill_value_for_type(type: pa.DataType) -> Any:
    if type not in DEFAULT_FILL_VALUE:
        raise ValueError("Unsupported type (do not know default fill)")

    return DEFAULT_FILL_VALUE[type]


def densendarray_datatype() -> st.SearchStrategy[pa.DataType]:
    # Arrow Tensor doesn't support bool_ or timestamp, and that is the only
    # read accessor we have. So for now, don't test those types.
    if HT_TEST_CONFIG["sc-61743_workaround"]:
        return ndarray_datatype().filter(lambda t: t != pa.bool_() and not pa.types.is_timestamp(t))

    return ndarray_datatype()


class SOMADenseNDArrayStateMachine(SOMANDArrayStateMachine):
    def __init__(self) -> None:
        super().__init__(shapes_factory=dense_array_shape)

    @initialize(type=densendarray_datatype(), shape=dense_array_shape())
    def setup(self, type, shape) -> None:
        super().setup(
            type,
            shape,
            soma.DenseNDArray.create(
                self.uri,
                type=type,
                shape=shape,
                context=self.context,
                tiledb_timestamp=None,  # TODO: no time-travel for now
            ),
        )

        # Initial state of dense ndarray should be completely filled with the
        # default TilDB fill value
        initial_array_state = pa.Tensor.from_numpy(
            np.full(
                self.shape,
                fill_value_for_type(self.type),
                dtype=self.type.to_pandas_dtype(),
            ),
        )
        assert initial_array_state.shape == self.shape
        assert initial_array_state.type == self.type

        # TODO: due to sc-61676, reads return incorrect results for any portion
        # of the array that has not be explicitly written.  Hack around by explicitly
        # writing fill values, AND disabling any resize operations.
        if HT_TEST_CONFIG["sc-61676_workaround"]:
            self.A.write(tuple(slice(0, n) for n in self.shape), initial_array_state)
            self._close()

        self.data_ledger = Ledger[ArrowTensorLedgerEntry](
            initial_entry=ArrowTensorLedgerEntry(
                data=initial_array_state,
                timestamp_ms=self.A.tiledb_timestamp_ms,
                name="initial entry",
            ),
            allows_duplicates=False,
        )

    def _array_exists(self, uri: str, context: soma.SOMAContext, tiledb_timestamp: int | None) -> bool:
        return soma.DenseNDArray.exists(uri, context=context, tiledb_timestamp=tiledb_timestamp)

    def _array_open(self, *, mode: OpenMode, tiledb_timestamp: int | None = None) -> None:
        self.A = soma.DenseNDArray.open(self.uri, mode=mode, context=self.context, tiledb_timestamp=tiledb_timestamp)

    ##
    # --- schema
    ##
    @precondition(lambda self: not self.closed)
    @invariant()
    def check_pytypes(self) -> None:
        assert isinstance(self.A, soma.DenseNDArray)
        assert self.A.soma_type == "SOMADenseNDArray"
        assert not self.A.is_sparse

    # XXX temporarily override this so we can disable any reshapes (sc-61676).
    # If we allow reshapes, read/write tests fail due to the bug.
    # TODO: remove this code (let base class do its thing) when this bug is fixed.
    @precondition(lambda self: HT_TEST_CONFIG["sc-61676_workaround"])
    @rule(data=st.data())
    def expand_shape(self, data: st.DataObject) -> None:
        return

    ##
    # --- data
    ##
    @precondition(lambda self: not self.closed and self.mode == "r")
    @rule(result_order=st.sampled_from(["row-major", "column-major"]))
    def check_read_all(self, result_order: str) -> None:
        tensor = self.A.read(result_order=result_order)
        expected = self.data_ledger.read(timestamp_ms=self.A.tiledb_timestamp_ms).to_tensor()
        if result_order != "row-major":
            expected = pa.Tensor.from_numpy(expected.to_numpy().T)

        assert self.type == tensor.type == expected.type
        assert tensor.shape == expected.shape
        if result_order != "row-major":
            assert tuple(reversed(tensor.shape)) == self.shape
        else:
            assert tensor.shape == self.shape
        assert np.array_equal(tensor.to_numpy(), expected.to_numpy(), equal_nan=True)

    @precondition(lambda self: not self.closed and self.mode == "r")
    @rule(data=st.data())
    def check_read_indexed(self, data: st.DataObject) -> None:
        inclusive_shape = tuple(s - 1 for s in self.shape)
        coords = data.draw(dense_indices(inclusive_shape))
        tensor = self.A.read(coords=coords).to_numpy()
        assert self.type.to_pandas_dtype() == tensor.dtype

        subslc = tuple(
            (slice(c.start, c.stop + 1 if c.stop is not None else None) if isinstance(c, slice) else slice(c, c + 1))
            for c in coords
        )
        expected = self.data_ledger.read(timestamp_ms=self.A.tiledb_timestamp_ms).to_numpy()[subslc]
        assert tensor.shape == expected.shape
        assert tensor.dtype == expected.dtype
        assert np.array_equal(tensor, expected, equal_nan=True)

    @precondition(lambda self: not self.closed and self.mode == "w")
    @precondition(
        lambda self: self.A.tiledb_timestamp_ms not in self.data_ledger.timestamps,
    )  # only one write per timestamp until sc-61223 and sc-61226 are fixed
    @rule(data=st.data())
    def write(self, data: st.DataObject) -> None:
        # draw sub-array
        ndim = len(self.shape)
        first = tuple(data.draw(st.integers(min_value=0, max_value=s - 1)) for s in self.shape)
        second = tuple(data.draw(st.integers(min_value=0, max_value=s - 1)) for s in self.shape)
        top_left = tuple(min(first[i], second[i]) for i in range(ndim))
        bot_right = tuple(max(first[i], second[i]) for i in range(ndim))
        coords = tuple(slice(top_left[i], bot_right[i]) for i in range(ndim))
        subarray = data.draw(
            ht_np.arrays(
                self.type.to_pandas_dtype(),
                shape=tuple(bot_right[i] - top_left[i] + 1 for i in range(ndim)),
            ),
        )

        # Write sub-array to the SOMA array
        fragments_before_write = get_entries(f"{self.uri}/__fragments")
        self.A.write(coords, pa.Tensor.from_numpy(subarray))
        new_fragments = set(get_entries(f"{self.uri}/__fragments")) - set(fragments_before_write)
        assert len(new_fragments) == 1

        # Save write in the ledger. The tensor ledger expects the "entire" array value,
        # not a differential value (i.e., it currently does not consolidate).
        merged_array = self.data_ledger.read(timestamp_ms=self.A.tiledb_timestamp_ms).to_numpy()
        inclusive_coords = tuple(slice(s.start, s.stop + 1) if isinstance(s, slice) else s for s in coords)
        merged_array[inclusive_coords] = subarray
        self.data_ledger.write(
            ArrowTensorLedgerEntry(
                timestamp_ms=self.A.tiledb_timestamp_ms,
                name=new_fragments.pop(),
                data=pa.Tensor.from_numpy(merged_array),
            ),
        )


TestSOMADenseNDArray = pytest.mark.usefixtures("setup_fixtures")(SOMADenseNDArrayStateMachine.TestCase)
