import datetime
import time
from unittest import mock

import pytest

import tiledbsoma
import tiledbsoma.options._soma_tiledb_context as stc
import tiledb


@pytest.fixture(autouse=True)
def global_ctx_reset():
    stc._default_global_ctx.cache_clear()
    yield


def test_lazy_init():
    """Verifies we don't construct a Ctx until we have to."""
    with mock.patch.object(tiledb, "Ctx", wraps=tiledb.Ctx) as mock_ctx:
        context = stc.SOMATileDBContext(tiledb_config={})
        assert context.tiledb_config == {
            "sm.mem.reader.sparse_global_order.ratio_array_data": 0.3
        }
        mock_ctx.assert_not_called()
        assert context._tiledb_ctx is None
        # Invoke the @property twice to ensure we only build one Ctx.
        with pytest.deprecated_call():
            assert context.tiledb_ctx is context.tiledb_ctx
        mock_ctx.assert_called_once()


def test_tiledb_ctx_init():
    config = {"hither": "yon"}
    with pytest.deprecated_call():
        context = stc.SOMATileDBContext(tiledb_ctx=tiledb.Ctx(config))
    assert "hither" in context.tiledb_config


def test_lazy_replace_config():
    """Verifies we don't construct a Ctx even if we call ``.replace``."""
    with mock.patch.object(tiledb, "Ctx", wraps=tiledb.Ctx) as mock_ctx:
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
    with pytest.deprecated_call():
        assert ctx.tiledb_ctx is ctx_2.tiledb_ctx


def test_unshared_ctx():
    """Verifies that contexts are not shared when not appropriate."""
    ctx = stc.SOMATileDBContext()
    ctx_2 = stc.SOMATileDBContext(tiledb_config={})
    ctx_3 = stc.SOMATileDBContext(tiledb_config={})
    with pytest.deprecated_call():
        assert ctx.tiledb_ctx is not ctx_2.tiledb_ctx
        assert ctx_2.tiledb_ctx is not ctx_3.tiledb_ctx


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


def test_replace_context():
    with pytest.deprecated_call():
        orig_ctx = stc.SOMATileDBContext(tiledb_ctx=tiledb.Ctx())
    new_tdb_ctx = tiledb.Ctx({"vfs.s3.region": "hy-central-1"})
    with pytest.deprecated_call():
        new_ctx = orig_ctx.replace(tiledb_ctx=new_tdb_ctx)
    with pytest.deprecated_call():
        assert new_ctx.tiledb_ctx is new_tdb_ctx


def test_replace_config_after_construction():
    context = stc.SOMATileDBContext()

    # verify defaults expected by subsequent tests
    assert context.timestamp_ms is None
    assert context.native_context.config()["vfs.s3.region"] == "us-east-1"

    now = int(time.time() * 1000)
    open_ts = context._open_timestamp_ms(None)
    assert -100 < now - open_ts < 100
    assert 999 == context._open_timestamp_ms(999)

    context_ts_1 = context.replace(timestamp=1)

    assert context_ts_1.timestamp_ms == 1
    assert context_ts_1._open_timestamp_ms(None) == 1
    assert context_ts_1._open_timestamp_ms(2) == 2

    with mock.patch.object(tiledb, "Ctx", wraps=tiledb.Ctx) as mock_ctx:
        # verify that the new context is lazily initialized.
        new_soma_ctx = context.replace(tiledb_config={"vfs.s3.region": "us-west-2"})
        assert new_soma_ctx.tiledb_config["vfs.s3.region"] == "us-west-2"
        mock_ctx.assert_not_called()
        with pytest.deprecated_call():
            new_tdb_ctx = new_soma_ctx.tiledb_ctx
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
