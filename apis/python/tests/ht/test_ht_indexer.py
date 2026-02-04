"""Hypothesis tests for IntIndexer module."""

from typing import Any, Union

import hypothesis as ht
import hypothesis.extra.numpy as ht_np
import numpy as np
import numpy.typing as npt
import pyarrow as pa
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

import tiledbsoma as soma
from tiledbsoma import pytiledbsoma as clib

from tests.ht._ht_util import (
    arrow_array_fast,
    arrow_chunked_array_fast,
    everything_except,
)


@given(
    data=ht_np.arrays(
        dtype=np.int64,
        shape=ht_np.array_shapes(max_dims=1, max_side=127),
        unique=True,
    ),
    context=st.one_of(st.from_type(soma.SOMAContext), st.none()),
)
def test_IntIndexer_ndarray_lookup(data: npt.NDArray[Any], context: soma.SOMAContext) -> None:
    assert np.array_equal(
        soma.IntIndexer(data=data, context=context).get_indexer(data),
        np.arange(0, len(data), dtype=np.int64),
    )


@given(
    data=st.one_of(
        (
            arrow_array_fast(np.int64, shape=st.integers(min_value=0, max_value=2047), unique=True),
            arrow_chunked_array_fast(
                dtype=np.int64,
                shape=st.integers(min_value=0, max_value=1023),
                splits=3,
                unique=True,
            ),
        ),
    ),
)
@settings(suppress_health_check=(ht.HealthCheck.function_scoped_fixture,))
def test_IntIndexer_arrow_lookup(data: pa.ChunkedArray, context: soma.SOMAContext) -> None:
    assert np.array_equal(
        soma.IntIndexer(data=data, context=context).get_indexer(data),
        np.arange(0, len(data), dtype=np.int64),
    )


@given(data=st.from_type(Union[np.ndarray[Any, Any], list[int]]))
@settings(suppress_health_check=(ht.HealthCheck.function_scoped_fixture,))
def test_fuzz_IntIndexer(data: npt.NDArray[Any], context: soma.SOMAContext) -> None:
    if isinstance(data, list):
        ht.assume(len(data) > 0 and any(not isinstance(x, int) for x in data))
    elif isinstance(data, np.ndarray):
        ht.assume(not (data.ndim == 1 and data.dtype == np.int64))
    with pytest.raises(Exception):
        soma.IntIndexer(data=data, context=context)


@given(
    data=ht_np.arrays(
        dtype=np.int64,
        shape=ht_np.array_shapes(max_dims=1, max_side=127),
        unique=True,
    ),
)
@settings(suppress_health_check=(ht.HealthCheck.function_scoped_fixture,))
def test_pytiledbsoma_IntIndexer_map_locations(data: npt.NDArray[np.int64], context: soma.SOMAContext) -> None:
    indexer = clib.IntIndexer(context._handle)
    indexer.map_locations(data)


@given(
    data=st.one_of(
        (
            ht_np.arrays(dtype=ht_np.array_dtypes(), shape=ht_np.array_shapes(), unique=True),
            ht_np.arrays(dtype=ht_np.array_dtypes(), shape=ht_np.array_shapes(), unique=False),
            st.from_type(Union[float, list, dict, str, bytearray]),
        ),
    ),
)
@settings(suppress_health_check=(ht.HealthCheck.function_scoped_fixture,))
def test_fuzz_pytiledbsoma_IntIndexer_map_locations(data: npt.NDArray[Any], context: soma.SOMAContext) -> None:
    ht.assume((not isinstance(data, np.ndarray)) or data.dtype != np.int64 or data.ndim != 1)

    indexer = clib.IntIndexer(context._handle)
    with pytest.raises(Exception):
        indexer.map_locations(data)


@given(
    data=st.one_of(
        (
            ht_np.arrays(dtype=ht_np.array_dtypes(), shape=ht_np.array_shapes(), unique=False),
            everything_except(np.ndarray),
        ),
    ),
)
@settings(suppress_health_check=(ht.HealthCheck.function_scoped_fixture,))
def test_fuzz_pytiledbsoma_Indexer_get_indexer_general(data: Any, context: soma.SOMAContext) -> None:
    ht.assume((not isinstance(data, np.ndarray)) or data.dtype != np.int64 or data.ndim != 1)

    indexer = clib.IntIndexer(context._handle)
    indexer.map_locations(np.arange(0, 100, dtype=np.int64))
    with pytest.raises(Exception):
        indexer.get_indexer_general(data)
