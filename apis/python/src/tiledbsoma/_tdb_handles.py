# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""Abstractions to more easily manage read and write access to TileDB data."""

from __future__ import annotations

from collections.abc import Iterator, MutableMapping, Sequence
from typing import Any, TypeVar, Union

import attrs
import numpy as np
from typing_extensions import Self

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


@attrs.define(frozen=True)
class MetadataWrapper(MutableMapping[str, Any]):
    """A wrapper storing the metadata of some TileDB object.

    Because the view of metadata does not change after open time, we immediately
    cache all of it and use that to handle all reads. Writes are then proxied
    through to the backing store and the cache is updated to match.
    """

    _owner: RawHandle
    """Tracks the modifications we have made to cache entries."""

    @classmethod
    def from_handle(cls, owner: RawHandle) -> Self:
        return cls(owner)

    def __contains__(self, key: object) -> bool:
        if self._owner.closed:
            raise SOMAError("Cannot get length of metadata; object is closed. Open in mode='r'.")
        return bool(self._owner.has_metadata(key))

    def __len__(self) -> int:
        if self._owner.closed:
            raise SOMAError("Cannot get length of metadata; object is closed. Open in mode='r'.")
        return int(self._owner.metadata_num())

    def __iter__(self) -> Iterator[str]:
        if self._owner.closed:
            raise SOMAError("Cannot iterate over metadata; object is closed. Open in mode='r'.")
        return iter(self._owner.meta)  # TODO: Add binding to expose the iterator

    def __getitem__(self, key: str) -> Any:  # noqa: ANN401
        if self._owner.closed:
            raise SOMAError(f"Cannot get metadata item '{key}'; object is closed. Open in mode='r'.")

        if not self._owner.has_metadata(key):
            raise KeyError

        return self._owner.get_metadata(key)

    def __setitem__(self, key: str, value: Any) -> None:  # noqa: ANN401
        if self._owner.closed:
            raise SOMAError(f"Cannot set metadata item '{key}'='{value}'; object is closed.")
        if self._owner.mode != "w":
            raise SOMAError(
                f"Cannot write metadata item to '{key}'; current mode='{self._owner.mode}'. Reopen in mode='w'."
            )
        _check_metadata_type(key, value)
        if isinstance(value, str):
            self._owner.set_metadata(key, np.array([value.encode("UTF-8")], "S"))
        elif isinstance(value, bytes):
            self._owner.set_metadata(key, np.array([value], "V"))
        else:
            self._owner.set_metadata(key, np.array([value]))

    def __delitem__(self, key: str) -> None:
        if self._owner.closed:
            raise SOMAError(f"Cannot delete metadata item at '{key}'; object is closed.")
        if self._owner.mode != "w":
            raise SOMAError(
                f"Cannot delete metadata item at '{key}'; current mode='{self._owner.mode}'. Reopen in mode='w'."
            )
        self._owner.delete_metadata(key)

    def __repr__(self) -> str:
        if self._owner.closed:
            return f"<{type(self).__name__}(<{type(self._owner).__name__} (closed)>)>"
        return f"<{type(self).__name__}(<{type(self._owner).__name__} (self.owner.mode)>) {self._owner.meta}>"


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
