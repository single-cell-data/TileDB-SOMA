# Copyright (c) 2021-2024 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2024 TileDB, Inc.
#
# Licensed under the MIT License.

"""Abstractions to more easily manage read and write access to TileDB data.

``open``, ``ArrayWrapper.open``, ``GroupWrapper.open`` are the important parts.
"""

import abc
import enum
from typing import (
    Any,
    Dict,
    Generic,
    Iterator,
    Mapping,
    MutableMapping,
    Optional,
    Tuple,
    Type,
    TypeVar,
    Union,
    cast,
)

import attrs
import numpy as np
from somacore import options
from typing_extensions import Literal, Self

from . import pytiledbsoma as clib
from ._exception import DoesNotExistError, SOMAError, is_does_not_exist_error
from ._types import METADATA_TYPES, Metadatum, OpenTimestamp
from .options._soma_tiledb_context import SOMATileDBContext

class _DictMod(enum.Enum):
    """State machine to keep track of modifications to a dictionary.

    This whole thing is a hack to allow users to treat the metadata dict
    like an actual dictionary because tiledb currently does not support multiple
    modifications to the same key (e.g., add-then-delete a metadata entry has
    undesired results) [sc-25089].
    """

    # Initially-absent keys are either added or not (added then removed).
    ABSENT = enum.auto()
    """The key is not present in the dict. Initial state."""
    ADDED = enum.auto()
    """The key was originally ABSENT but has been added."""

    # Initially-present keys can be either updated or deleted.
    PRESENT = enum.auto()
    """The key is in the dict and is unchanged. Initial state."""
    UPDATED = enum.auto()
    """The key was originally PRESENT but has been changed."""
    DELETED = enum.auto()
    """The key was originally PRESENT but has been deleted."""

    @classmethod
    def start_state(cls, dct: Mapping[Any, Any], key: Any) -> "_DictMod":
        """Returns the starting state for a DictMod given the key of dct."""
        return cls.PRESENT if key in dct else cls.ABSENT

    def next_state(self, action: Literal["set", "del"]) -> "_DictMod":
        """Determines the next state of an entry given the action."""
        return {
            _DictMod.ABSENT: {
                "set": _DictMod.ADDED,
            },
            _DictMod.ADDED: {
                "set": _DictMod.ADDED,
                "del": _DictMod.ABSENT,
            },
            _DictMod.PRESENT: {
                "set": _DictMod.UPDATED,
                "del": _DictMod.DELETED,
            },
            _DictMod.UPDATED: {
                "set": _DictMod.UPDATED,
                "del": _DictMod.DELETED,
            },
            _DictMod.DELETED: {
                "set": _DictMod.UPDATED,
            },
        }[self][action]


@attrs.define(frozen=True)
class MetadataWrapper(MutableMapping[str, Any]):
    """A wrapper storing the metadata of some TileDB object.

    Because the view of metadata does not change after open time, we immediately
    cache all of it and use that to handle all reads. Writes are then proxied
    through to the backing store and the cache is updated to match.
    """

    # XXX TEMP
    # owner: Wrapper[CLibHandle]
    owner: Any

    cache: Dict[str, Any]
    _mods: Dict[str, "_DictMod"] = attrs.field(init=False, factory=dict)
    """Tracks the modifications we have made to cache entries."""

    def __len__(self) -> int:
        self.owner._check_open()
        return len(self.cache)

    def __iter__(self) -> Iterator[str]:
        self.owner._check_open()
        return iter(self.cache)

    def __getitem__(self, key: str) -> Any:
        self.owner._check_open()
        return self.cache[key]

    def __setitem__(self, key: str, value: Any) -> None:
        self.owner.writer  # Ensures we're open in write mode.
        state = self._current_state(key)
        _check_metadata_type(key, value)
        self.cache[key] = value
        self._mods[key] = state.next_state("set")

    def __delitem__(self, key: str) -> None:
        self.owner.writer  # Ensures we're open in write mode.
        state = self._current_state(key)
        del self.cache[key]
        self._mods[key] = state.next_state("del")

    def _current_state(self, key: str) -> _DictMod:
        return self._mods.get(key, _DictMod.start_state(self.cache, key))

    def _write(self) -> None:
        """Writes out metadata changes, if there were any."""
        if not self._mods:
            # There were no changes (e.g., it's a read handle).  Do nothing.
            return
        # Only try to get the writer if there are changes to be made.
        for key, mod in self._mods.items():
            if mod in (_DictMod.ADDED, _DictMod.UPDATED):
                set_metadata = self.owner._handle.set_metadata
                val = self.cache[key]
                if isinstance(val, str):
                    set_metadata(key, np.array([val], "S"))
                else:
                    set_metadata(key, np.array([val]))
            if mod is _DictMod.DELETED:
                self.owner._handle.delete_metadata(key)

        # Temporary hack: When we flush writes, note that the cache
        # is back in sync with disk.
        self._mods.clear()

    def __repr__(self) -> str:
        prefix = f"{type(self).__name__}({self.owner})"
        if self.owner.closed:
            return f"<{prefix}>"
        return f"<{prefix} {self.cache}>"


def _check_metadata_type(key: str, obj: Metadatum) -> None:
    """Pre-checks that a metadata entry can be stored in an array.

    These checks are reproduced from the TileDB Python metadata-setting methods,
    but are slightly more restrictive than what TileDB allows in general:
    TileDB allows (some) arrays as metadata values, but the SOMA spec does not
    allow arrays of any kind.

    We have to pre-check since we don't write metadata changes until closing.
    """
    if not isinstance(key, str):
        raise TypeError(f"metadata keys must be strings, not {type(key)}")
    if isinstance(obj, METADATA_TYPES):
        return
    raise TypeError(f"cannot store {type(obj)} instance as metadata")
