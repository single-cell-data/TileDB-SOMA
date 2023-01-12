from typing import Union

import pytest
import tiledb

from tiledbsoma import TileDBObject, TileDBSessionContext


class TestTDBObject(TileDBObject):
    @property
    def soma_type(self) -> str:
        return self.__class__.__name__

    def _tiledb_open(self, mode: str = "r") -> Union[tiledb.Array, tiledb.Group]:
        return tiledb.Group("")


def test_tiledb_object_default_session_context_added():
    o = TestTDBObject("parent")
    assert o._soma_session_context is not None
    assert o._soma_session_context.tiledb_ctx is not None
    assert o._ctx is not None
    assert o._ctx == o._soma_session_context.tiledb_ctx


def test_child_inherits_parent_session_context():
    session_context = TileDBSessionContext()
    p = TestTDBObject("parent", session_context=session_context)

    c = TestTDBObject("parent/child", parent=p)

    assert c._soma_session_context is not None
    assert c._soma_session_context == p._soma_session_context


def test_mutually_exclusive_create_args():
    session_context = TileDBSessionContext()
    p = TestTDBObject("parent", session_context=session_context)

    with pytest.raises(AssertionError):
        TestTDBObject("parent/child", parent=p, session_context=session_context)


