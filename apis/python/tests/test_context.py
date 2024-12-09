import datetime
import importlib
import time
from unittest import mock

import pytest

import tiledbsoma
import tiledbsoma.options._soma_tiledb_context as stc
import tiledbsoma.pytiledbsoma as clib


def test_lazy_init():
    """Verifies we don't construct a Ctx until we have to."""
    with mock.patch.object(clib, "SOMAContext", wraps=clib.SOMAContext) as mock_ctx:
        context = stc.SOMATileDBContext(tiledb_config={})
        assert context.tiledb_config == {
            "sm.mem.reader.sparse_global_order.ratio_array_data": 0.3
        }
        mock_ctx.assert_not_called()
        assert context._native_context is None
        # Invoke the @property twice to ensure we only build one Ctx.
        assert context.native_context is not None
        assert context.native_context is context.native_context
        mock_ctx.assert_called_once()


def test_lazy_replace_config():
    """Verifies we don't construct a Ctx even if we call ``.replace``."""
    with mock.patch.object(clib, "SOMAContext", wraps=clib.SOMAContext) as mock_ctx:
        context = stc.SOMATileDBContext()
        new_context = context.replace(tiledb_config={"hello": "goodbye"})
        assert new_context.tiledb_config == {
            "sm.mem.reader.sparse_global_order.ratio_array_data": 0.3,
            "hello": "goodbye",
        }
        mock_ctx.assert_not_called()


def test_delete_config_entry():
    context = stc.SOMATileDBContext(tiledb_config={"hither": "yon"})
    assert context.tiledb_config == {
        "hither": "yon",
        "sm.mem.reader.sparse_global_order.ratio_array_data": 0.3,
    }
    new_context = context.replace(tiledb_config={"hither": None})
    # We've removed the only non-default entry; this should work.
    assert new_context.tiledb_config == {
        "sm.mem.reader.sparse_global_order.ratio_array_data": 0.3
    }


def test_shared_ctx():
    """Verifies that one global context is shared by default."""
    ctx = stc.SOMATileDBContext()
    ctx_2 = stc.SOMATileDBContext()
    assert ctx.native_context is ctx_2.native_context


def test_unshared_ctx():
    """Verifies that contexts are not shared when not appropriate."""
    ctx = stc.SOMATileDBContext()
    ctx_2 = stc.SOMATileDBContext(tiledb_config={})
    ctx_3 = stc.SOMATileDBContext(tiledb_config={})
    assert ctx.native_context is not ctx_2.native_context
    assert ctx_2.native_context is not ctx_3.native_context


def test_replace_timestamp():
    orig_ctx = stc.SOMATileDBContext()
    assert orig_ctx.timestamp is None
    assert orig_ctx.timestamp_ms is None
    ts_ctx = orig_ctx.replace(timestamp=1683817200000)
    assert ts_ctx.timestamp == datetime.datetime(
        2023, 5, 11, 15, 0, tzinfo=datetime.timezone.utc
    )
    assert ts_ctx.timestamp_ms == 1683817200000
    same_ts_ctx = ts_ctx.replace()  # replace nothing!
    assert ts_ctx.timestamp == same_ts_ctx.timestamp
    no_ts_ctx = ts_ctx.replace(timestamp=None)
    assert no_ts_ctx.timestamp is None


def test_replace_config_after_construction():
    context = stc.SOMATileDBContext()

    # verify defaults expected by subsequent tests
    assert context.timestamp_ms is None
    if tiledbsoma.pytiledbsoma.embedded_version_triple() < (2, 27, 0):
        assert context.native_context.config()["vfs.s3.region"] == "us-east-1"
    else:
        assert context.native_context.config()["vfs.s3.region"] == ""

    now = int(time.time() * 1000)
    open_ts = context._open_timestamp_ms(None)
    assert -100 < now - open_ts < 100
    assert 999 == context._open_timestamp_ms(999)

    context_ts_1 = context.replace(timestamp=1)

    assert context_ts_1.timestamp_ms == 1
    assert context_ts_1._open_timestamp_ms(None) == 1
    assert context_ts_1._open_timestamp_ms(2) == 2

    with mock.patch.object(clib, "SOMAContext", wraps=clib.SOMAContext) as mock_ctx:
        # verify that the new context is lazily initialized.
        new_soma_ctx = context.replace(tiledb_config={"vfs.s3.region": "us-west-2"})
        assert new_soma_ctx.tiledb_config["vfs.s3.region"] == "us-west-2"
        mock_ctx.assert_not_called()
        new_tdb_ctx = new_soma_ctx.native_context
        mock_ctx.assert_called_once()
        assert new_tdb_ctx.config()["vfs.s3.region"] == "us-west-2"


def test_malformed_concurrency_config_value():
    import numpy as np

    with pytest.raises(tiledbsoma.SOMAError):
        ctx = tiledbsoma.SOMATileDBContext(
            tiledb_config={"soma.compute_concurrency_level": "not-a-number"}
        )

        tiledbsoma.IntIndexer(np.arange(100, dtype=np.int64), context=ctx).get_indexer(
            np.array([0, 1])
        )


@pytest.mark.skipif(
    importlib.util.find_spec("tiledb") is not None, reason="TileDB-Py is installed"
)
def test_tiledb_ctx_without_tiledb():
    # Test that tiledb_ctx errors out as expected without tiledb-py

    with pytest.raises(ModuleNotFoundError):
        tiledbsoma.SOMATileDBContext(tiledb_ctx="junk")

    sctx = tiledbsoma.SOMATileDBContext()
    with pytest.raises(ModuleNotFoundError):
        sctx.tiledb_ctx

    with pytest.raises(ModuleNotFoundError):
        sctx.replace(tiledb_ctx="junk")


@pytest.mark.skipif(
    importlib.util.find_spec("tiledb") is None, reason="TileDB-Py is not installed"
)
def test_tiledb_ctx_with_tiledb():
    # If tiledb-py is installed, test that tiledb_ctx works to handle tiledb.Ctx
    import tiledb

    # Default
    sctx = tiledbsoma.SOMATileDBContext(tiledb_ctx=tiledb.Ctx())
    assert sctx.tiledb_ctx.config() == tiledb.Ctx().config()

    # Pass config
    sctx = tiledbsoma.SOMATileDBContext(tiledb_ctx=tiledb.Ctx({"foo": "bar"}))
    assert sctx.tiledb_ctx.config()["foo"] == "bar"

    # Replace config
    sctx = sctx.replace(tiledb_ctx=tiledb.Ctx({"foo": "baz"}))
    assert sctx.tiledb_ctx.config()["foo"] == "baz"
