from typing import Union

import tiledb

from tiledbsoma import TileDBObject


class FakeTileDBObject(TileDBObject):
    mock_open_error: bool = False

    @property
    def soma_type(self) -> str:
        return self.__class__.__name__

    @property
    def _tiledb_obj(self) -> Union[tiledb.Array, tiledb.Group]:
        return tiledb.Group("")

    def _sub_open(self) -> None:
        if self.mock_open_error:
            raise RuntimeError("mock")


def test_tiledb_object_default_context_added():
    o = FakeTileDBObject("parent")
    assert o.context is not None
    assert o.context.tiledb_ctx is not None
    assert o._ctx is not None
    assert o._ctx == o.context.tiledb_ctx


def test_open_close():
    o = FakeTileDBObject("parent")

    assert not o.mode
    o.open_legacy("r")
    assert o.mode == "r"
    o.close()
    assert not o.mode
    o.close()  # no-op
    assert not o.mode

    with o.open_legacy("w"):
        assert o.mode == "w"
    assert not o.mode

    try:
        with o.open_legacy("w"):
            assert o.mode == "w"
            raise RuntimeError("test")
    except Exception:
        pass
    assert not o.mode

    # Must open to use as contextmanager
    bad = False
    try:
        with o:
            bad = True
        bad = True
    except Exception:
        pass
    assert not bad

    # Internal _ensure_open()
    with o._ensure_open():
        assert o.mode == "r"
    assert not o.mode
    with o.open_legacy():
        with o._ensure_open():
            assert o.mode == "r"
        assert o.mode == "r"
    assert not o.mode
    bad = False
    try:
        # reject attempt to write when already open read-only
        with o.open_legacy("r"):
            with o._ensure_open("w"):
                bad = True
            bad = True
    except Exception:
        pass
    assert not bad
    assert not o.mode
    with o.open_legacy("w"):
        # allow attempt to read when already open for write (manly for metadata; underlying TileDB
        # objects may reject other attempts)
        with o._ensure_open("r"):
            assert o.mode == "w"
        assert o.mode == "w"
    assert not o.mode

    # closed after an exception during open()
    bad = False
    o.mock_open_error = True
    try:
        with o.open_legacy("r"):
            bad = True
        bad = True
    except Exception:
        pass
    assert not bad
    assert o.closed
