from __future__ import annotations

import numpy as np
import pytest

from tiledbsoma import IntIndexer, SOMAContext


@pytest.mark.parametrize("context", [None, SOMAContext()])
def test_reindexer_api(context: SOMAContext | None):
    keys = np.arange(3, 10, 2)
    ids = np.arange(3, 10, 2)
    expected = np.array([0, 1, 2, 3])
    indexer = IntIndexer(keys, context=context)
    result = indexer.get_indexer(ids)
    assert np.equal(result.all(), expected.all())
