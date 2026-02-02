# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.
"""Type and interface declarations that are not specific to options.

These declarations and functions are not directly part of the user-facing API,
but are intended to be used by SOMA implementations for annotations and
their own internal type-checking purposes.
"""

from __future__ import annotations

from collections.abc import Sequence
from concurrent import futures
from typing import TypeVar

from typing_extensions import Protocol, TypeGuard

StatusAndReason = tuple[bool, str]
"""Information for whether an upgrade-shape or resize would succeed
if attempted, along with a reason why not."""


def is_nonstringy_sequence(it: object) -> bool:
    """Returns true if a sequence is a "normal" sequence and not str or bytes.

    str and bytes are "weird" sequences because iterating them gives you
    another str or bytes instance for each character, and when used as a
    sequence is not what users want.
    """
    return not isinstance(it, (str, bytes)) and isinstance(it, Sequence)


def to_string_tuple(obj: str | Sequence[str]) -> tuple[str, ...]:
    """Returns a tuple of string values.

    If the input is a string, it is returned as a tuple with the string as its
    only item. If it is otherwise a sequence of strings, the sequence is converted
    to a tuple.
    """
    return (obj,) if isinstance(obj, str) else tuple(obj)


def str_or_seq_length(obj: str | Sequence[str]) -> int:
    """Returns the number of str values.

    If input is a string, returns 1. Otherwise, returns the number of strings in the
    sequence.
    """
    return 1 if isinstance(obj, str) else len(obj)


_T = TypeVar("_T")
_T_co = TypeVar("_T_co", covariant=True)


class Slice(Protocol[_T_co]):
    """A slice which stores a certain type of object.

    This protocol describes the built-in ``slice`` type, with a hint to callers
    about what type they should put *inside* the slice.  It is for type
    annotations only and is not runtime-checkable (i.e., you can't do
    ``isinstance(thing, Slice)``), because ``range`` objects also have
    ``start``/``stop``/``step`` and would match, but are *not* slices.
    """

    @property
    def start(self) -> _T_co | None: ...

    @property
    def stop(self) -> _T_co | None: ...

    @property
    def step(self) -> _T_co | None: ...


def is_slice_of(__obj: object, __typ: type[_T]) -> TypeGuard[Slice[_T]]:
    return (
        # We only respect `slice`s proper.
        isinstance(__obj, slice)
        and (__obj.start is None or isinstance(__obj.start, __typ))
        and (__obj.stop is None or isinstance(__obj.stop, __typ))
        and (__obj.step is None or isinstance(__obj.step, __typ))
    )


class ContextBase(Protocol):
    """A protocol for a context manager that can be used as a base class.
    If a threadpool is specified as part of the context, it will be used by
    experiment queries. Otherwise, the implementer will use its own threadpool.
    """

    threadpool: futures.ThreadPoolExecutor | None
