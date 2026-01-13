from unittest import mock

import pytest

import tiledbsoma
import tiledbsoma.pytiledbsoma as clib


def test_lazy_init():
    """Verifies we don't construct a context until we have to."""
    with mock.patch.object(clib, "SOMAContext", wraps=clib.SOMAContext) as mock_ctx:
        context = tiledbsoma.SOMAContext(config={})
        assert context.config == {"sm.mem.reader.sparse_global_order.ratio_array_data": "0.3"}
        mock_ctx.assert_not_called()
        assert context._native_context is None
        # Invoke the @property twice to ensure we only build one context.
        assert context.native_context is not None
        assert context.native_context is context.native_context
        mock_ctx.assert_called_once()


def test_lazy_replace_config():
    """Verifies we don't construct a context even if we call ``.replace``."""
    with mock.patch.object(clib, "SOMAContext", wraps=clib.SOMAContext) as mock_ctx:
        context = tiledbsoma.SOMAContext()
        new_context = context.replace(config={"hello": "goodbye"})
        assert new_context.config == {
            "sm.mem.reader.sparse_global_order.ratio_array_data": "0.3",
            "hello": "goodbye",
        }
        mock_ctx.assert_not_called()


def test_delete_config_entry():
    context = tiledbsoma.SOMAContext(
        config={
            "hither": "yon",
            "sm.mem.reader.sparse_global_order.ratio_array_data": "0.5",
        }
    )
    assert context.config == {
        "hither": "yon",
        "sm.mem.reader.sparse_global_order.ratio_array_data": "0.5",
    }
    new_context = context.replace(config={"hither": None})
    assert new_context.config == {"sm.mem.reader.sparse_global_order.ratio_array_data": "0.5"}


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
    assert ctx1.native_context is ctx2.native_context
    assert ctx2.native_context is ctx3.native_context


def test_default_context_with_config():
    """Verifies that contexts are not shared when directly constructed."""
    tiledbsoma.SOMAContext._default_context = None
    ctx1 = tiledbsoma.SOMAContext.set_default(config={"hither": "yon"})
    ctx2 = tiledbsoma.SOMAContext.get_default()
    ctx3 = tiledbsoma.SOMAContext.get_default()
    assert ctx1.native_context is ctx2.native_context
    assert ctx2.native_context is ctx3.native_context


def test_change_default_context():
    tiledbsoma.SOMAContext._default_context = None
    ctx1 = tiledbsoma.SOMAContext.set_default()
    ctx2 = tiledbsoma.SOMAContext.set_default(replace=True)
    assert ctx1.native_context is not ctx2.native_context


def test_not_shared_ctx():
    """Verifies that contexts are not shared when not appropriate."""
    tiledbsoma.SOMAContext._default_context = None
    tiledbsoma.SOMAContext.set_default()
    ctx1 = tiledbsoma.SOMAContext.get_default()
    ctx2 = tiledbsoma.SOMAContext()
    ctx3 = tiledbsoma.SOMAContext()
    assert ctx1.native_context is not ctx2.native_context
    assert ctx2.native_context is not ctx3.native_context
    assert ctx1.native_context is not ctx3.native_context


def test_replace_config_after_construction():
    context = tiledbsoma.SOMAContext()

    # verify defaults expected by subsequent tests
    if tiledbsoma.pytiledbsoma.tiledb_version() < (2, 27, 0):
        assert context.native_context.config()["vfs.s3.region"] == "us-east-1"
    else:
        assert not context.native_context.config()["vfs.s3.region"]

    with mock.patch.object(clib, "SOMAContext", wraps=clib.SOMAContext) as mock_ctx:
        # verify that the new context is lazily initialized.
        new_soma_ctx = context.replace(config={"vfs.s3.region": "us-west-2"})
        assert new_soma_ctx.config["vfs.s3.region"] == "us-west-2"
        mock_ctx.assert_not_called()
        new_tdb_ctx = new_soma_ctx.native_context
        mock_ctx.assert_called_once()
        assert new_tdb_ctx.config()["vfs.s3.region"] == "us-west-2"


def test_malformed_concurrency_config_value():
    import numpy as np

    with pytest.raises(tiledbsoma.SOMAError):
        ctx = tiledbsoma.SOMAContext(config={"soma.compute_concurrency_level": "not-a-number"})

        tiledbsoma.IntIndexer(np.arange(100, dtype=np.int64), context=ctx).get_indexer(np.array([0, 1]))
