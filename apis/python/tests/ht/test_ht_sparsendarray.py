"""Hypothesis tests for SparseNDArray."""

from __future__ import annotations

import datetime
import shutil
from collections.abc import Mapping, Sequence
from typing import Any, Union

import hypothesis as ht
import numpy as np
import pyarrow as pa
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from hypothesis.stateful import (
    initialize,
    invariant,
    precondition,
    rule,
)
from somacore.options import OpenMode

import tiledbsoma as soma
import tiledbsoma._sparse_nd_array

from tests.ht._array_state_machine import SOMANDArrayStateMachine
from tests.ht._ht_util import (
    arrow_array,
    ndarray_datatype,
    posix_filename,
    tables_equal,
    tiledb_timestamps,
)
from tests.ht._ledger import ArrowTableLedgerEntry, Ledger, get_entries


@st.composite
def sparse_array_shape(
    draw: st.DrawFn,
    *,
    min_shape: tuple[int, ...] | None = None,
    max_shape: tuple[int, ...] | None = None,
    allow_none: bool = False,
) -> tuple[int | None, ...]:
    """Strategy to generate nd array shapes."""

    MAX_DIMS = 7
    if min_shape is not None:
        ndim = len(min_shape)
    elif max_shape is not None:
        ndim = len(max_shape)
    else:
        ndim = draw(st.integers(min_value=1, max_value=MAX_DIMS))

    min_values = [1] * ndim if min_shape is None else min_shape

    # due to how we draw random coordinates in sparse_array() strategy, the product
    # of the shape must fit in an int64. SOMA/TileDB has some weird internal restrictions
    # that force the max of any one dimension to be a bit smaller. Set our per-dim
    # limit to the min of the nth root of int64.max or the tiledb limit.
    shape_limit = min(int((2**63 - 1) ** (1 / ndim)), (2**63 - 2050))
    max_values = [shape_limit] * ndim if max_shape is None else [min(shape_limit, s) for s in max_shape]

    if allow_none:
        elements = [
            draw(
                st.one_of(
                    st.none(),
                    st.integers(min_value=min_values[i], max_value=max_values[i]),
                ),
            )
            for i in range(ndim)
        ]
    else:
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
        assert np.prod(new_shape) <= 2**63 - 1

    return new_shape


@st.composite
def sparse_array(
    draw: st.DrawFn,
    shape: tuple[int, ...],
    schema: pa.Schema,
    *,
    density: float | None = None,
) -> pa.Table:
    """
    Draw a sparse array with ndim, in SOMA Table format, with types matching
    schema, and dimensions within the given shape.
    """
    MAX_NNZ = 1 * 1024**2  # caps memory use

    shape_prod = np.prod(shape)
    if shape_prod == 0:
        return schema.empty_table()

    if density is None:
        max_density = MAX_NNZ / shape_prod if shape_prod > MAX_NNZ else 1.0
        assert 0 <= max_density <= 1
        density = draw(st.floats(min_value=0, max_value=max_density))

    nnz = max(1, int(density * shape_prod))
    assert nnz <= MAX_NNZ, "Sparse array is too large."

    rng = np.random.default_rng(seed=draw(st.integers(min_value=0)))
    coords = rng.choice(shape_prod, size=nnz, replace=False)
    tbl_dict = {}
    for n, n_len in reversed(list(enumerate(shape))):
        coords, tbl_dict[f"soma_dim_{n}"] = np.divmod(coords, n_len)
    assert np.all(coords == 0)
    tbl_dict = dict(reversed(tbl_dict.items()))

    type = schema.field("soma_data").type
    tbl_dict["soma_data"] = draw(arrow_array(type, shape=nnz))
    return pa.Table.from_pydict(tbl_dict, schema=schema)


@given(
    uri=posix_filename(),
    type=st.from_type(pa.DataType).filter(
        lambda t: (
            pa.types.is_primitive(t)
            and not (pa.types.is_timestamp(t) and t.tz is not None)
            and not pa.types.is_time(t)
            and not pa.types.is_date(t)
            and t
            not in [
                pa.float16(),
            ]
        ),
    ),
    shape=st.lists(
        st.one_of(st.none(), st.integers(min_value=1, max_value=2**31 - 1)),
        min_size=1,
        max_size=10,
    ),
    platform_config=st.from_type(Union[dict[str, str], None]),
    context=st.from_type(Union[soma.SOMAContext, None]),
    tiledb_timestamp=tiledb_timestamps(),
)
@settings(suppress_health_check=(ht.HealthCheck.function_scoped_fixture,))
def test_fuzz_SparseNDArray_create(
    tmp_path,
    uri: str,
    type: pa.DataType,
    shape: Sequence[int | None],
    platform_config: dict[str, Mapping[str, Any]] | object | None,
    context: tiledbsoma.SOMAContext | None,
    tiledb_timestamp: int | datetime.datetime | None,
) -> None:
    try:
        fname = (tmp_path / uri).as_posix()
        if None in shape:
            with pytest.deprecated_call():
                A = soma.SparseNDArray.create(
                    uri=fname,
                    type=type,
                    shape=shape,
                    platform_config=platform_config,
                    context=context,
                    tiledb_timestamp=tiledb_timestamp,
                )
        else:
            A = soma.SparseNDArray.create(
                uri=fname,
                type=type,
                shape=shape,
                platform_config=platform_config,
                context=context,
                tiledb_timestamp=tiledb_timestamp,
            )
        A.close()

        with soma.open(fname, context=context) as A:
            assert len(A.schema.types) == len(shape) + 1
            assert A.schema.field("soma_data").type == type
            assert A.shape == tuple((s or 1) for s in shape)
            assert A.soma_type == "SOMASparseNDArray"

    finally:
        shutil.rmtree(tmp_path / uri, ignore_errors=True)


class SOMASparseNDArrayStateMachine(SOMANDArrayStateMachine):
    def __init__(self) -> None:
        super().__init__(shapes_factory=sparse_array_shape)

    @initialize(type=ndarray_datatype(), shape=sparse_array_shape(allow_none=False))
    def setup(self, type: pa.DataType, shape: tuple[int | None, ...]) -> None:
        if None in shape:
            with pytest.deprecated_call():
                super().setup(
                    type,
                    shape,
                    soma.SparseNDArray.create(
                        self.uri,
                        type=type,
                        shape=shape,
                        context=self.context,
                        tiledb_timestamp=None,  # no time-travel for now
                    ),
                )
        else:
            super().setup(
                type,
                shape,
                soma.SparseNDArray.create(
                    self.uri,
                    type=type,
                    shape=shape,
                    context=self.context,
                    tiledb_timestamp=None,  # no time-travel for now
                ),
            )

        self.data_ledger = Ledger[ArrowTableLedgerEntry](
            initial_entry=ArrowTableLedgerEntry(
                data=self.schema.empty_table(),
                timestamp_ms=self.A.tiledb_timestamp_ms,
                name="initial entry",
                index_columns=[f"soma_dim_{n}" for n in range(len(shape))],
            ),
            allows_duplicates=False,
        )

    def _array_exists(self, uri: str, context: soma.SOMAContext, tiledb_timestamp: int | None) -> bool:
        return soma.SparseNDArray.exists(uri, context=context, tiledb_timestamp=tiledb_timestamp)

    def _array_open(self, *, mode: OpenMode, tiledb_timestamp: int | None = None) -> None:
        self.A = soma.SparseNDArray.open(self.uri, mode=mode, context=self.context, tiledb_timestamp=tiledb_timestamp)

    ##
    # --- schema
    ##
    @precondition(lambda self: not self.closed)
    @invariant()
    def check_pytypes(self) -> None:
        assert isinstance(self.A, soma.SparseNDArray)
        assert self.A.soma_type == "SOMASparseNDArray"
        assert self.A.is_sparse

    ##
    # --- data
    ##

    @precondition(lambda self: not self.closed and self.mode == "r")
    @invariant()
    def check_read_all(self) -> None:
        timestamp_ms = self.A.tiledb_timestamp_ms
        sort_order = [(f"soma_dim_{n}", "ascending") for n in range(len(self.shape))]
        expected = self.data_ledger.read(timestamp_ms=timestamp_ms).to_table().sort_by(sort_order)
        found = self.A.read().tables().concat().sort_by(sort_order)
        assert tables_equal(
            found,
            expected,
            equal_nan=bool(pa.types.is_floating(self.type)),
        ), f"{found}\n is not equal to {expected}"

    @precondition(lambda self: not self.closed and self.mode == "r")
    @invariant()
    def check_nnz(self) -> None:
        expected = len(self.data_ledger.read(timestamp_ms=self.A.tiledb_timestamp_ms).to_table())
        assert expected == self.A.nnz, "NNZ mismatch"
        assert self.A.nnz <= self.A._handle.fragment_cell_count(), "NNZ vs fragment_cell_count inconsistency"

    @precondition(lambda self: not self.closed and self.mode == "w")
    @precondition(
        lambda self: self.A.tiledb_timestamp_ms not in self.data_ledger.timestamps,
    )  # only one write per timestamp until sc-61223 and sc-61226 are fixed
    @rule(data=st.data())
    def write(self, data: st.DataObject) -> None:
        coo_tbl = data.draw(sparse_array(self.shape, self.schema))

        fragments_before_write = get_entries(f"{self.uri}/__fragments")
        self.A.write(coo_tbl)
        new_fragments = set(get_entries(f"{self.uri}/__fragments")) - set(fragments_before_write)
        assert len(new_fragments) == 1
        self.data_ledger.write(
            ArrowTableLedgerEntry(
                timestamp_ms=self.A.tiledb_timestamp_ms,
                name=new_fragments.pop(),
                data=coo_tbl,
                index_columns=[f"soma_dim_{n}" for n in range(len(self.shape))],
            ),
        )


TestSOMASparseNDArray = pytest.mark.usefixtures("setup_fixtures")(SOMASparseNDArrayStateMachine.TestCase)
