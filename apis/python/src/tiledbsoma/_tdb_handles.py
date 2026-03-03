# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""Abstractions to more easily manage read and write access to TileDB data."""

from __future__ import annotations

import enum
from collections.abc import Iterator, Mapping, MutableMapping, Sequence
from typing import Any, TypeVar, Union

import attrs
import numpy as np
from typing_extensions import Literal, Self

from . import pytiledbsoma as clib
from ._core_options import OpenMode
from ._exception import SOMAError
from ._types import METADATA_TYPES, Metadatum

AxisDomain = Union[tuple[Any, Any], list[Any], None]
Domain = Sequence[AxisDomain]

RawHandle = Union[
    clib.SOMAArray,
    clib.SOMADataFrame,
    clib.SOMAPointCloudDataFrame,
    clib.SOMAGeometryDataFrame,
    clib.SOMASparseNDArray,
    clib.SOMADenseNDArray,
    clib.SOMAGroup,
    clib.SOMACollection,
    clib.SOMAMeasurement,
    clib.SOMAExperiment,
    clib.SOMAScene,
    clib.SOMAMultiscaleImage,
]
_RawHdl_co = TypeVar("_RawHdl_co", bound=RawHandle, covariant=True)
"""A raw TileDB object. Covariant because Handles are immutable enough."""

_SOMAObjectType = TypeVar("_SOMAObjectType", bound=clib.SOMAObject)


def _open_mode_to_clib_mode(mode: OpenMode) -> clib.OpenMode:
    """Convert OpenMode to clib.OpenMode."""
    if mode == "r":
        return clib.OpenMode.soma_read
    if mode == "w":
        return clib.OpenMode.soma_write
    if mode == "d":
        return clib.OpenMode.soma_delete
    raise ValueError(f"Unexpected mode '{mode}'. Valid modes are 'r', 'w', or 'd'.")


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
    def start_state(cls, dct: Mapping[Any, Any], key: Any) -> _DictMod:  # noqa: ANN401
        """Returns the starting state for a DictMod given the key of dct."""
        return cls.PRESENT if key in dct else cls.ABSENT

    def next_state(self, action: Literal["set", "del"]) -> _DictMod:
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

    owner: RawHandle
    cache: dict[str, Any]
    _mods: dict[str, _DictMod] = attrs.field(init=False, factory=dict)
    """Tracks the modifications we have made to cache entries."""

    @classmethod
    def from_handle(cls, owner: RawHandle) -> Self:
        return cls(owner, dict(owner.meta))

    def __len__(self) -> int:
        if self.owner.closed:
            raise SOMAError("Cannot get length of metadata; object is closed. Open in mode='r'.")
        return len(self.cache)

    def __iter__(self) -> Iterator[str]:
        if self.owner.closed:
            raise SOMAError("Cannot iterate over metadata; object is closed. Open in mode='r'.")
        return iter(self.cache)

    def __getitem__(self, key: str) -> Any:  # noqa: ANN401
        if self.owner.closed:
            raise SOMAError(f"Cannot get metadata item '{key}'; object is closed. Open in mode='r'.")
        return self.cache[key]

    def __setitem__(self, key: str, value: Any) -> None:  # noqa: ANN401
        if self.owner.closed:
            raise SOMAError(f"Cannot set metadata item '{key}'='{value}'; object is closed.")
        if self.owner.mode != "w":
            raise SOMAError(
                f"Cannot write metadata item to '{key}'; current mode='{self.owner.mode}'. Reopen in mode='w'."
            )
        state = self._current_state(key)
        _check_metadata_type(key, value)
        self.cache[key] = value
        self._mods[key] = state.next_state("set")

    def __delitem__(self, key: str) -> None:
        if self.owner.closed:
            raise SOMAError(f"Cannot delete metadata item at '{key}'; object is closed.")
        if self.owner.mode != "w":
            raise SOMAError(
                f"Cannot delete metadata item at '{key}'; current mode='{self.owner.mode}'. Reopen in mode='w'."
            )
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

        errors: list[Exception] = []
        for key, mod in self._mods.items():
            try:
                if mod in (_DictMod.ADDED, _DictMod.UPDATED):
                    val = self.cache[key]
                    if isinstance(val, str):
                        self.owner.set_metadata(key, np.array([val.encode("UTF-8")], "S"))
                    elif isinstance(val, bytes):
                        self.owner.set_metadata(key, np.array([val], "V"))
                    else:
                        self.owner.set_metadata(key, np.array([val]))
                if mod is _DictMod.DELETED:
                    self.owner.delete_metadata(key)
            except Exception as e:  # noqa: BLE001, PERF203
                # This should be done with Exception Groups
                errors.append(e)

        # Temporary hack: When we flush writes, note that the cache
        # is back in sync with disk.
        self._mods.clear()

        if errors:
            details = [repr(error) for error in errors]

            error_msg_details = "\n".join(details)

            raise SOMAError(
                f"[MetadataWrapper][_write] {len(errors)} error(s) occured while writing metadata to disk. Details: \n {error_msg_details}",
            )

    def __repr__(self) -> str:
        if self.owner.closed:
            return f"<{type(self).__name__}(<{type(self.owner).__name__} (closed)>)>"
        return f"<{type(self).__name__}(<{type(self.owner).__name__} (self.owner.mode)>) {self.cache}>"


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
