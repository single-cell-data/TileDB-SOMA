from __future__ import annotations

from typing import Any

import hypothesis as ht
import hypothesis.extra.numpy as ht_np
import hypothesis.strategies as st
import pyarrow as pa
import pytest

import tiledbsoma as soma

from tests.ht._ht_test_config import HT_TEST_CONFIG
from tests.ht._ht_util import (
    arrow_array_fast,
    arrow_chunked_array_fast,
    arrow_datatypes,
    arrow_shape,
)


@pytest.fixture(scope="class")
def setup_fixtures(request, tmp_path_factory, soma_tiledb_context: soma.SOMAContext) -> None:
    """Set a class variable pointing at required fixtures - useful for Hypothesis RuleBasedStateMachine test objects."""
    request.cls.tmp_path_factory = tmp_path_factory
    request.cls.soma_tiledb_context = soma_tiledb_context


@pytest.fixture
def ht_test_config() -> dict[str, Any]:
    return HT_TEST_CONFIG


@pytest.fixture
def context() -> soma.SOMAContext:
    return soma.SOMAContext()


# Register Hypothesis strategies for use with `strategies.from_type()`
st.register_type_strategy(pa.DataType, arrow_datatypes())
st.register_type_strategy(
    pa.Array,
    arrow_array_fast(
        dtype=ht_np.array_dtypes(),
        shape=arrow_shape(shape=st.integers(min_value=0, max_value=2047)),
    ),
)
st.register_type_strategy(
    pa.ChunkedArray,
    arrow_chunked_array_fast(
        dtype=ht_np.array_dtypes(),
        shape=arrow_shape(shape=st.integers(min_value=0, max_value=4095)),
    ),
)
# TODO: vary context configuration?
st.register_type_strategy(soma.SOMAContext, st.just(soma.SOMAContext()))


# Register hypothesis profile for extensive/expensive test runs
ht.settings.register_profile(
    "expensive",
    max_examples=10000,
    print_blob=True,
)
