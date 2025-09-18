# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""Abstractions to more easily manage read and write access to TileDB data.

``open``, ``ArrayWrapper.open``, ``GroupWrapper.open`` are the important parts.
"""

from __future__ import annotations

import abc
import enum
import warnings
from collections.abc import Iterator, Mapping, MutableMapping, Sequence
from typing import (
    Any,
    Generic,
    TypeVar,
    Union,
    cast,
)

import attrs
import numpy as np
import pyarrow as pa
from somacore import options
from typing_extensions import Literal, Self

from . import pytiledbsoma as clib
from ._constants import (
    SOMA_ENCODING_VERSION_METADATA_KEY,
    SOMA_OBJECT_TYPE_METADATA_KEY,
    SUPPORTED_SOMA_ENCODING_VERSIONS,
)
from ._exception import DoesNotExistError, SOMAError, is_does_not_exist_error
from ._types import METADATA_TYPES, Metadatum, OpenTimestamp
from .options._soma_tiledb_context import SOMATileDBContext

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


def _open_mode_to_clib_mode(mode: options.OpenMode) -> clib.OpenMode:
    """Convert options.OpenMode to clib.OpenMode."""
    if mode == "r":
        return clib.OpenMode.soma_read
    if mode == "w":
        return clib.OpenMode.soma_write
    if mode == "d":
        return clib.OpenMode.soma_delete
    raise ValueError(f"Unexpected mode '{mode}'. Valid modes are 'r', 'w', or 'd'.")


def open_handle_wrapper(
    uri: str,
    mode: options.OpenMode,
    context: SOMATileDBContext,
    timestamp: OpenTimestamp | None,
    clib_type: str | None = None,
) -> Wrapper[RawHandle]:
    """Determine whether the URI is an array or group, and open it."""
    timestamp_ms = context._open_timestamp_ms(timestamp)

    type_to_class_ = {
        "somadataframe": DataFrameWrapper,
        "somapointclouddataframe": PointCloudDataFrameWrapper,
        "somageometrydataframe": GeometryDataFrameWrapper,
        "somadensendarray": DenseNDArrayWrapper,
        "somasparsendarray": SparseNDArrayWrapper,
        "somacollection": CollectionWrapper,
        "somaexperiment": ExperimentWrapper,
        "somameasurement": MeasurementWrapper,
        "somascene": SceneWrapper,
        "somamultiscaleimage": MultiscaleImageWrapper,
    }

    if clib_type is None or clib_type.lower() in ["somaarray", "somagroup"]:
        try:
            open_mode = _open_mode_to_clib_mode(mode)
            handle = clib.SOMAObject.open(
                uri=uri,
                mode=open_mode,
                context=context.native_context,
                timestamp=(0, timestamp_ms),
                clib_type=clib_type,
            )
        except Exception as tdbe:
            if is_does_not_exist_error(tdbe):
                raise DoesNotExistError(tdbe) from tdbe
            raise
        try:
            return type_to_class_[handle.type.lower()].open_from_handle(handle, uri=uri, mode=mode, context=context)
        except KeyError:
            raise SOMAError(f"{uri!r} has unknown storage type {clib_type!r}") from None

    try:
        return type_to_class_[clib_type.lower()].open(uri=uri, mode=mode, context=context, timestamp=timestamp_ms)
    except KeyError:
        raise SOMAError(f"{uri!r} has unknown storage type {clib_type!r}") from None


@attrs.define(eq=False, hash=False, slots=False)
class Wrapper(Generic[_RawHdl_co], metaclass=abc.ABCMeta):
    """Wrapper for TileDB handles to manage lifecycle and metadata.

    Callers may read and use (non-underscored) members but should never set
    attributes on instances.
    """

    uri: str
    mode: options.OpenMode
    context: SOMATileDBContext
    timestamp_ms: int
    _handle: _RawHdl_co
    closed: bool = attrs.field(default=False, init=False)
    clib_type: str | None = None

    @classmethod
    def open(
        cls,
        uri: str,
        mode: options.OpenMode,
        context: SOMATileDBContext,
        timestamp: OpenTimestamp | None,
    ) -> Self:
        if mode not in ("r", "w", "d"):
            raise ValueError(f"Invalid open mode {mode!r}")
        timestamp_ms = context._open_timestamp_ms(timestamp)

        try:
            tdb_handle = cls._opener(uri, mode, context, timestamp_ms)
        except (RuntimeError, SOMAError) as tdbe:
            if is_does_not_exist_error(tdbe):
                raise DoesNotExistError(tdbe) from tdbe
            raise SOMAError(tdbe) from tdbe
        return cls.open_from_handle(tdb_handle, uri=uri, mode=mode, context=context)

    @classmethod
    def open_from_handle(
        cls, clib_handle: clib.SOMAObject, *, uri: str, mode: options.OpenMode, context: SOMATileDBContext
    ) -> Self:
        timestamp = context._open_timestamp_ms(clib_handle.timestamp)
        try:
            wrapper = cls(uri, mode, context, timestamp, clib_handle)
        except RuntimeError as tdbe:
            if is_does_not_exist_error(tdbe):
                raise DoesNotExistError(tdbe) from tdbe
            raise
        wrapper._do_initial_reads(clib_handle)
        obj_type = wrapper.metadata.get(SOMA_OBJECT_TYPE_METADATA_KEY)
        if obj_type is None:
            raise SOMAError(
                f"Cannot access stored TileDB object with TileDB-SOMA. The object is missing "
                f"the required '{SOMA_OBJECT_TYPE_METADATA_KEY!r}' metadata key.",
            )
        if isinstance(obj_type, bytes):
            obj_type = str(obj_type, "utf-8")
        if not isinstance(obj_type, str):
            raise SOMAError(
                f"Cannot access stored TileDB object with TileDB-SOMA. The metadata key "
                f"'{SOMA_OBJECT_TYPE_METADATA_KEY!r}' has unexpected type '{type(obj_type)}'.",
            )

        encoding_version = wrapper.metadata.get(SOMA_ENCODING_VERSION_METADATA_KEY)
        if encoding_version is None:
            raise SOMAError(
                f"Cannot access stored TileDB object with TileDB-SOMA. The object is missing "
                f"the required '{SOMA_ENCODING_VERSION_METADATA_KEY!r}' metadata key.",
            )
        if isinstance(encoding_version, bytes):
            encoding_version = str(encoding_version, "utf-8")
        if encoding_version not in SUPPORTED_SOMA_ENCODING_VERSIONS:
            raise ValueError(
                f"Unsupported SOMA object encoding version '{encoding_version}'. TileDB-SOMA "
                f"needs to be updated to a more recent version.",
            )

        return wrapper

    @classmethod
    @abc.abstractmethod
    def _opener(
        cls,
        uri: str,
        mode: options.OpenMode,
        context: SOMATileDBContext,
        timestamp: int,
    ) -> _RawHdl_co:
        """Opens and returns a TileDB object specific to this type."""
        raise NotImplementedError

    # Covariant types should normally not be in parameters, but this is for
    # internal use only so it's OK.
    def _do_initial_reads(self, reader: _RawHdl_co) -> None:  # type: ignore[misc]
        """Final setup step before returning the Handle.

        This is passed a raw TileDB object opened in read mode, since writers
        will need to retrieve data from the backing store on setup.
        """
        # non-attrs-managed field
        self.metadata = MetadataWrapper(self, dict(reader.meta))

    @property
    def reader(self) -> _RawHdl_co:
        """Accessor to assert that you are working in read mode."""
        if self.closed:
            raise SOMAError(f"{self} is closed")
        if self.mode == "r":
            return self._handle
        raise SOMAError(f"Cannot read from {self}; current mode='{self.mode}'. Reopen in mode='r'.")

    @property
    def writer(self) -> _RawHdl_co:
        """Accessor to assert that you are working in write mode."""
        if self.closed:
            raise SOMAError(f"{self} is closed")
        if self.mode == "w":
            return self._handle
        raise SOMAError(f"Cannot write to {self}; current mode='{self.mode}'. Reopen in mode='w'.")

    @property
    def deleter(self) -> _RawHdl_co:
        """Accessor to assert that you are working in delete mode."""
        if self.closed:
            raise SOMAError(f"{self} is closed")
        if self.mode == "d":
            return self._handle
        if self.mode == "w":
            warnings.warn(
                f"Deleting in write mode is deprecated. {self} should be reopened with mode='d'.",
                DeprecationWarning,
                stacklevel=3,
            )
            return self._handle
        raise SOMAError(f"Cannot delete from {self}; current mode='{self.mode}'. Reopen in mode='d'.")  # noqa: S608

    def close(self) -> None:
        if self.closed:
            return
        self.metadata._write()
        self._handle.close()
        self.closed = True

    def _check_open(self) -> None:
        if self.closed:
            raise SOMAError(f"{self!r} is closed")

    def __repr__(self) -> str:
        closed_str = " (closed)" if self.closed else ""
        return f"<{type(self).__name__} {self.mode} on {self.uri!r}{closed_str}>"

    def __enter__(self) -> Self:
        return self

    def __exit__(self, *_: Any) -> None:  # noqa: ANN401
        self.close()

    def __del__(self) -> None:
        self.close()


AnyWrapper = Wrapper[RawHandle]
"""Non-instantiable type representing any Handle."""


@attrs.define(frozen=True)
class GroupEntry:
    uri: str
    wrapper_type: type[AnyWrapper]

    @classmethod
    def from_soma_group_entry(cls, obj: tuple[str, str]) -> GroupEntry:
        uri, type = obj[0], obj[1]
        if type == "SOMAArray":
            return GroupEntry(uri, SOMAArrayWrapper)
        if type == "SOMAGroup":
            return GroupEntry(uri, SOMAGroupWrapper)
        raise SOMAError(f"internal error: unknown object type {uri}")


class SOMAGroupWrapper(Wrapper[_SOMAObjectType]):
    """Base class for Pybind11 SOMAGroupWrapper handles."""

    _WRAPPED_TYPE: type[_SOMAObjectType]

    clib_type = "SOMAGroup"

    @classmethod
    def _opener(
        cls,
        uri: str,
        mode: options.OpenMode,
        context: SOMATileDBContext,
        timestamp: int,
    ) -> clib.SOMAGroup:
        open_mode = _open_mode_to_clib_mode(mode)
        return cls._WRAPPED_TYPE.open(
            uri,
            mode=open_mode,
            context=context.native_context,
            timestamp=(0, timestamp),
        )

    def _do_initial_reads(self, group: clib.SOMAGroup) -> None:
        super()._do_initial_reads(group)

        self.initial_contents = {
            name: GroupEntry.from_soma_group_entry(entry) for name, entry in group.members().items()
        }

    @property
    def meta(self) -> MetadataWrapper:
        return self.metadata

    def members(self) -> dict[str, tuple[str, str]]:
        return cast("dict[str, tuple[str, str]]", self._handle.members())


class CollectionWrapper(SOMAGroupWrapper[clib.SOMACollection]):
    """Wrapper around a Pybind11 CollectionWrapper handle."""

    _WRAPPED_TYPE = clib.SOMACollection


class ExperimentWrapper(SOMAGroupWrapper[clib.SOMAExperiment]):
    """Wrapper around a Pybind11 ExperimentWrapper handle."""

    _WRAPPED_TYPE = clib.SOMAExperiment


class MeasurementWrapper(SOMAGroupWrapper[clib.SOMAMeasurement]):
    """Wrapper around a Pybind11 MeasurementWrapper handle."""

    _WRAPPED_TYPE = clib.SOMAMeasurement


class MultiscaleImageWrapper(SOMAGroupWrapper[clib.SOMAMultiscaleImage]):
    """Wrapper around a Pybind11 MultiscaleImage handle."""

    _WRAPPED_TYPE = clib.SOMAMultiscaleImage


class SceneWrapper(SOMAGroupWrapper[clib.SOMAScene]):
    """Wrapper around a Pybind11 SceneWrapper handle."""

    _WRAPPED_TYPE = clib.SOMAScene


class SOMAArrayWrapper(Wrapper[_SOMAObjectType]):
    """Base class for Pybind11 SOMAArrayWrapper handles."""

    _WRAPPED_TYPE: type[_SOMAObjectType]

    clib_type = "SOMAArray"

    @classmethod
    def _opener(
        cls,
        uri: str,
        mode: options.OpenMode,
        context: SOMATileDBContext,
        timestamp: int,
    ) -> clib.SOMAArray:
        open_mode = _open_mode_to_clib_mode(mode)
        return cls._WRAPPED_TYPE.open(
            uri,
            mode=open_mode,
            context=context.native_context,
            timestamp=(0, timestamp),
        )

    def _do_initial_reads(self, reader: RawHandle) -> None:
        """Final setup step before returning the Handle.

        This is passed a raw TileDB object opened in read mode, since writers
        will need to retrieve data from the backing store on setup.
        """
        # non-attrs-managed field
        self.metadata = MetadataWrapper(self, dict(reader.meta))

    @property
    def schema(self) -> pa.Schema:
        return self._handle.schema

    def schema_config_options(self) -> clib.PlatformSchemaConfig:
        """Returns a class containing the TileDB platform configuration options that
        can be read from an array schema.
        """
        return self._handle.schema_config_options()

    @property
    def meta(self) -> MetadataWrapper:
        return self.metadata

    @property
    def ndim(self) -> int:
        return len(self._handle.dimension_names)

    @property
    def shape(self) -> tuple[int, ...]:
        """Not implemented for DataFrame."""
        return cast("tuple[int, ...]", tuple(self._handle.shape))

    @property
    def maxshape(self) -> tuple[int, ...]:
        """Not implemented for DataFrame."""
        return cast("tuple[int, ...]", tuple(self._handle.maxshape))


class DataFrameWrapper(SOMAArrayWrapper[clib.SOMADataFrame]):
    """Wrapper around a Pybind11 SOMADataFrame handle."""

    _WRAPPED_TYPE = clib.SOMADataFrame


class PointCloudDataFrameWrapper(SOMAArrayWrapper[clib.SOMAPointCloudDataFrame]):
    """Wrapper around a Pybind11 SOMAPointCloudDataFrame handle."""

    _WRAPPED_TYPE = clib.SOMAPointCloudDataFrame


class GeometryDataFrameWrapper(SOMAArrayWrapper[clib.SOMAGeometryDataFrame]):
    """Wrapper around a Pybind11 SOMAGeometryDataFrame handle."""

    _WRAPPED_TYPE = clib.SOMAGeometryDataFrame


class DenseNDArrayWrapper(SOMAArrayWrapper[clib.SOMADenseNDArray]):
    """Wrapper around a Pybind11 DenseNDArrayWrapper handle."""

    _WRAPPED_TYPE = clib.SOMADenseNDArray


class SparseNDArrayWrapper(SOMAArrayWrapper[clib.SOMASparseNDArray]):
    """Wrapper around a Pybind11 SparseNDArrayWrapper handle."""

    _WRAPPED_TYPE = clib.SOMASparseNDArray


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

    owner: Wrapper[RawHandle]
    cache: dict[str, Any]
    _mods: dict[str, _DictMod] = attrs.field(init=False, factory=dict)
    """Tracks the modifications we have made to cache entries."""

    def __len__(self) -> int:
        self.owner._check_open()
        return len(self.cache)

    def __iter__(self) -> Iterator[str]:
        self.owner._check_open()
        return iter(self.cache)

    def __getitem__(self, key: str) -> Any:  # noqa: ANN401
        self.owner._check_open()
        return self.cache[key]

    def __setitem__(self, key: str, value: Any) -> None:  # noqa: ANN401
        self.owner.writer  # noqa: B018 Ensures we're open in write mode.
        state = self._current_state(key)
        _check_metadata_type(key, value)
        self.cache[key] = value
        self._mods[key] = state.next_state("set")

    def __delitem__(self, key: str) -> None:
        self.owner.writer  # noqa: B018 Ensures we're open in write mode.
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
                    set_metadata = self.owner._handle.set_metadata
                    val = self.cache[key]
                    if isinstance(val, str):
                        set_metadata(key, np.array([val.encode("UTF-8")], "S"))
                    elif isinstance(val, bytes):
                        set_metadata(key, np.array([val], "V"))
                    else:
                        set_metadata(key, np.array([val]))
                if mod is _DictMod.DELETED:
                    self.owner._handle.delete_metadata(key)
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
