from typing import Union

import pytest
import tiledb

from tiledbsoma import TileDBObject
from tiledbsoma.options import SOMATileDBContext


class FakeTileDBObject(TileDBObject):
    @property
    def soma_type(self) -> str:
        return self.__class__.__name__

    def _tiledb_open(self, mode: str = "r") -> Union[tiledb.Array, tiledb.Group]:
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

    with pytest.raises(AssertionError):
        FakeTileDBObject("parent/child", parent=p, context=context)
