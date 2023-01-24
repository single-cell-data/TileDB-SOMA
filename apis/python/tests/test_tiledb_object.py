from typing import Union

import pytest
import tiledb

from tiledbsoma import TileDBObject
from tiledbsoma.options import SOMATileDBContext


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


def test_child_inherits_parent_context():
    context = SOMATileDBContext()
    p = FakeTileDBObject("parent", context=context)

    c = FakeTileDBObject("parent/child", parent=p)

    assert c.context is not None
    assert c.context == p.context


def test_mutually_exclusive_create_args():
    context = SOMATileDBContext()
    p = FakeTileDBObject("parent", context=context)

    with pytest.raises(TypeError):
        FakeTileDBObject("child", parent=p, context=context)


def test_open_close():
    o = FakeTileDBObject("parent")

    assert not o.mode
    o.open("r")
    assert o.mode == "r"
    o.close()
    assert not o.mode
    o.close()  # no-op
    assert not o.mode

    with o.open("w"):
        assert o.mode == "w"
    assert not o.mode

    try:
        with o.open("w"):
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
    with o.open():
        with o._ensure_open():
            assert o.mode == "r"
        assert o.mode == "r"
    assert not o.mode
    bad = False
    try:
        # reject attempt to write when already open read-only
        with o.open("r"):
            with o._ensure_open("w"):
                bad = True
            bad = True
    except Exception:
        pass
    assert not bad
    assert not o.mode
    with o.open("w"):
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
        with o.open("r"):
            bad = True
        bad = True
    except Exception:
        pass
    assert not bad
    assert o.closed
