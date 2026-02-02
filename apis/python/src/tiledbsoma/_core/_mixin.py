# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.
"""Tools for making mixins with SOMA Collections."""

from __future__ import annotations

from collections.abc import MutableMapping
from typing import Generic, TypeVar, overload

import attrs

from . import base

_ST = TypeVar("_ST", bound=base.SOMAObject)
_Coll = MutableMapping[str, _ST]
_T = TypeVar("_T")


@attrs.define()
class item(Generic[_T]):
    """Descriptor to transform property access into indexing.

    This descriptor works on mapping objects to allow simple specification of
    properties that are backed by map entries::

        class FirstSecondMixin:

            first = item(str)
            second = item(int, "2nd")

        class FSCollection(FirstSecondMixin, CollectionBase):
            pass

        inst = FSCollection(...)

        # This is equivalent to getting inst["first"]
        inst.first

        # This is equivalent to setting inst["2nd"]
        inst.second = 500
    """

    typ: type[_T] | None = None
    """The type we expect to return from this field."""

    item_name: str | None = None
    """The name of the item we are getting (``x._backing["whatever"]``).

    This uses the name of the field by default but can be manually overridden.
    """

    field_name: str = attrs.field(default="<unknown>", init=False)
    """The name of this field (``x.whatever``). Set automatically."""

    def __set_name__(self, owner: type[_Coll], name: str) -> None:  # type: ignore[type-arg]
        del owner  # unused
        self.field_name = name
        if self.item_name is None:
            self.item_name = name

    @overload
    def __get__(self, inst: None, owner: type[_Coll]) -> item[_T]: ...  # type: ignore[type-arg]

    @overload
    def __get__(self, inst: _Coll, owner: type[_Coll]) -> _T: ...  # type: ignore[type-arg]

    def __get__(self, inst: _Coll | None, owner: type[_Coll]) -> item[_T] | _T:  # type: ignore[type-arg]
        del owner  # unused
        if not inst:
            return self
        assert self.item_name is not None
        try:
            # TODO: Type-check params/returns?
            return inst[self.item_name]  # type: ignore[no-any-return]
        except KeyError as ke:
            raise AttributeError(f"{_typename(inst)!r} object has no attribute {self.field_name!r}") from ke

    def __set__(self, inst: _Coll, value: _T) -> None:  # type: ignore[type-arg]
        assert self.item_name is not None
        # Pretend it's a MutableMapping for the type-checker.
        # If it fails that's OK; we need to raise anyway.
        try:
            inst[self.item_name] = value
        except KeyError as ke:
            raise AttributeError(f"{_typename(inst)!r} does not support assigning to item {self.item_name!r}") from ke

    def __delete__(self, inst: _Coll) -> None:  # type: ignore[type-arg]
        assert self.item_name is not None
        try:
            del inst[self.item_name]
        except KeyError as ke:
            raise AttributeError(f"{_typename(inst)} does not support deleting {self.item_name!r}") from ke


def _typename(x: object) -> str:
    return type(x).__name__
