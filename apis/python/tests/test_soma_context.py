import pytest

import tiledbsoma


def test_delete_config_entry():
    context = tiledbsoma.SOMAContext.create(
        config={
            "hither": "yon",
            "sm.mem.reader.sparse_global_order.ratio_array_data": "0.5",
        }
    )
    assert context.config["hither"] == "yon"
    assert context.config["sm.mem.reader.sparse_global_order.ratio_array_data"] == "0.5"
    new_context = context.replace(config={"hither": None})
    assert "hither" not in new_context.config
    assert new_context.config["sm.mem.reader.sparse_global_order.ratio_array_data"] == "0.5"


def test_default_context_not_set():
    tiledbsoma.SOMAContext._default_context = None
    with pytest.raises(RuntimeError):
        tiledbsoma.SOMAContext.get_default()


def test_default_conext_already_set():
    tiledbsoma.SOMAContext._default_context = None
    tiledbsoma.SOMAContext.set_default()
    with pytest.raises(RuntimeError):
        tiledbsoma.SOMAContext.set_default()


def test_default_context_no_config():
    """Verifies that contexts are not shared when directly constructed."""
    tiledbsoma.SOMAContext._default_context = None
    ctx1 = tiledbsoma.SOMAContext.set_default()
    ctx2 = tiledbsoma.SOMAContext.get_default()
    ctx3 = tiledbsoma.SOMAContext.get_default()
    assert ctx1._handle is ctx2._handle
    assert ctx2._handle is ctx3._handle


def test_default_context_with_config():
    """Verifies that contexts are not shared when directly constructed."""
    tiledbsoma.SOMAContext._default_context = None
    ctx1 = tiledbsoma.SOMAContext.set_default(config={"hither": "yon"})
    ctx2 = tiledbsoma.SOMAContext.get_default()
    ctx3 = tiledbsoma.SOMAContext.get_default()
    assert ctx1._handle is ctx2._handle
    assert ctx2._handle is ctx3._handle


def test_change_default_context():
    tiledbsoma.SOMAContext._default_context = None
    ctx1 = tiledbsoma.SOMAContext.set_default()
    ctx2 = tiledbsoma.SOMAContext.set_default(replace=True)
    assert ctx1._handle is not ctx2._handle


def test_not_shared_ctx():
    """Verifies that contexts are not shared when not appropriate."""
    tiledbsoma.SOMAContext._default_context = None
    tiledbsoma.SOMAContext.set_default()
    ctx1 = tiledbsoma.SOMAContext.get_default()
    ctx2 = tiledbsoma.SOMAContext.create()
    ctx3 = tiledbsoma.SOMAContext.create()
    assert ctx1._handle is not ctx2._handle
    assert ctx2._handle is not ctx3._handle
    assert ctx1._handle is not ctx3._handle


def test_malformed_concurrency_config_value():
    import numpy as np

    with pytest.raises(tiledbsoma.SOMAError):
        ctx = tiledbsoma.SOMAContext.create(config={"soma.compute_concurrency_level": "not-a-number"})

        tiledbsoma.IntIndexer(np.arange(100, dtype=np.int64), context=ctx).get_indexer(np.array([0, 1]))
