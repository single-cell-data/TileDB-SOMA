from typing import Union

import pytest
import tiledb

from tiledbsoma import TileDBObject
from tiledbsoma.options import SOMATileDBContext


class FakeTileDBObject(TileDBObject):
    @property
    def soma_type(self) -> str:
        return self.__class__.__name__

    @property
    def _tiledb_obj(self) -> Union[tiledb.Array, tiledb.Group]:
        return tiledb.Group("")


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

    assert not o._open_mode
    o.open("r")
    assert o._open_mode == "r"
    o.close()
    assert not o._open_mode
    o.close()  # no-op
    assert not o._open_mode

    with o.open("w"):
        assert o._open_mode == "w"
    assert not o._open_mode

    try:
        with o.open("w"):
            assert o._open_mode == "w"
            raise RuntimeError("test")
        assert False
    except:  # noqa: E722
        assert not o._open_mode

    # Must open to use as contextmanager
    bad = False
    try:
        with o:
            bad = True
        bad = True
    except:  # noqa: E722
        pass
    assert not bad
