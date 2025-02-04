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
def make_tmp_dir(request, tmp_path_factory) -> None:
    """Set a class variable - useful for Hypothesis RuleBasedStateMachine test objects."""
    request.cls.tmp_path_factory = tmp_path_factory


@pytest.fixture
def ht_test_config() -> dict[str, Any]:
    return HT_TEST_CONFIG


@pytest.fixture
def context() -> soma.SOMATileDBContext:
    return soma.SOMATileDBContext()


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
st.register_type_strategy(soma.SOMATileDBContext, st.just(soma.SOMATileDBContext()))


# Register hypothesis profile for extensive/expensive test runs
ht.settings.register_profile("expensive", max_examples=10000, print_blob=True)
