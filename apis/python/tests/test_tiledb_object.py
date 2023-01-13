from typing import Union

import pytest
import tiledb

from tiledbsoma import TileDBObject, SomaTileDBContext


class TestTDBObject(TileDBObject):
    @property
    def soma_type(self) -> str:
        return self.__class__.__name__

    def _tiledb_open(self, mode: str = "r") -> Union[tiledb.Array, tiledb.Group]:
        return tiledb.Group("")


def test_tiledb_object_default_context_added():
    o = TestTDBObject("parent")
    assert o.context is not None
    assert o.context.tiledb_ctx is not None
    assert o._ctx is not None
    assert o._ctx == o.context.tiledb_ctx


def test_child_inherits_parent_context():
    context = SomaTileDBContext()
    p = TestTDBObject("parent", context=context)

    c = TestTDBObject("parent/child", parent=p)

    assert c.context is not None
    assert c.context == p.context


def test_mutually_exclusive_create_args():
    context = SomaTileDBContext()
    p = TestTDBObject("parent", context=context)

    with pytest.raises(AssertionError):
        TestTDBObject("parent/child", parent=p, context=context)
