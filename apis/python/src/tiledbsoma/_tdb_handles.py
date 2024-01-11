# Copyright (c) 2021-2023 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2023 TileDB, Inc.
#
# Licensed under the MIT License.

"""Abstractions to more easily manage read and write access to TileDB data.

``open``, ``ArrayWrapper.open``, ``GroupWrapper.open`` are the important parts.
"""

import abc
import enum
from typing import (
    Any,
    Callable,
    Dict,
    Generic,
    Iterator,
    List,
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
import pyarrow as pa
import tiledb
from numpy.typing import DTypeLike
from somacore import options
from typing_extensions import Literal, Self

from . import pytiledbsoma as clib
from ._exception import DoesNotExistError, SOMAError, is_does_not_exist_error
from ._types import OpenTimestamp
from .options._soma_tiledb_context import SOMATileDBContext

RawHandle = Union[tiledb.Array, tiledb.Group, clib.SOMADataFrame]
_RawHdl_co = TypeVar("_RawHdl_co", bound=RawHandle, covariant=True)
"""A raw TileDB object. Covariant because Handles are immutable enough."""


def open(
    uri: str,
    mode: options.OpenMode,
    context: SOMATileDBContext,
    timestamp: Optional[OpenTimestamp],
) -> "Wrapper[RawHandle]":
    """Determine whether the URI is an array or group, and open it."""
    obj_type = tiledb.object_type(uri, ctx=context.tiledb_ctx)
    if not obj_type:
        raise DoesNotExistError(f"{uri!r} does not exist")

    try:
        return _open_with_clib_wrapper(uri, mode, context, timestamp)
    except SOMAError:
        # This object still uses tiledb-py and must be handled below
        if obj_type == "array":
            return ArrayWrapper.open(uri, mode, context, timestamp)
        if obj_type == "group":
            return GroupWrapper.open(uri, mode, context, timestamp)

    # Invalid object
    raise SOMAError(f"{uri!r} has unknown storage type {obj_type!r}")


def _open_with_clib_wrapper(
    uri: str,
    mode: options.OpenMode,
    context: SOMATileDBContext,
    timestamp: Optional[OpenTimestamp] = None,
) -> "DataFrameWrapper":
    open_mode = clib.OpenMode.read if mode == "r" else clib.OpenMode.write
    config = {k: str(v) for k, v in context.tiledb_config.items()}
    timestamp_ms = context._open_timestamp_ms(timestamp)
    obj = clib.SOMAObject.open(uri, open_mode, config, (0, timestamp_ms))
    if obj.type == "SOMADataFrame":
        return DataFrameWrapper._from_soma_object(obj, context)
    raise SOMAError(f"clib.SOMAObject {obj.type!r} not yet supported")


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

    @classmethod
    def open(
        cls,
        uri: str,
        mode: options.OpenMode,
        context: SOMATileDBContext,
        timestamp: Optional[OpenTimestamp],
    ) -> Self:
        if mode not in ("r", "w"):
            raise ValueError(f"Invalid open mode {mode!r}")
        timestamp_ms = context._open_timestamp_ms(timestamp)
        try:
            tdb = cls._opener(uri, mode, context, timestamp_ms)
            handle = cls(uri, mode, context, timestamp_ms, tdb)
            if mode == "w":
                with cls._opener(uri, "r", context, timestamp_ms) as auxiliary_reader:
                    handle._do_initial_reads(auxiliary_reader)
            else:
                handle._do_initial_reads(tdb)
        except tiledb.TileDBError as tdbe:
            if is_does_not_exist_error(tdbe):
                raise DoesNotExistError(f"{uri!r} does not exist") from tdbe
            raise
        return handle

    @classmethod
    def _from_soma_object(
        cls, soma_object: clib.SOMAObject, context: SOMATileDBContext
    ) -> Self:
        uri = soma_object.uri
        mode = soma_object.mode
        timestamp = soma_object.timestamp
        try:
            handle = cls(uri, mode, context, timestamp, soma_object)
            if handle.mode == "w":
                with cls._opener(uri, mode, context, timestamp) as auxiliary_reader:
                    handle._do_initial_reads(auxiliary_reader)
            else:
                handle._do_initial_reads(soma_object)
        except tiledb.TileDBError as tdbe:
            if is_does_not_exist_error(tdbe):
                raise DoesNotExistError(f"{handle.uri!r} does not exist") from tdbe
            raise
        return handle

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
        raise NotImplementedError()

    # Covariant types should normally not be in parameters, but this is for
    # internal use only so it's OK.
    def _do_initial_reads(self, reader: _RawHdl_co) -> None:  # type: ignore[misc]
        """Final setup step before returning the Handle.

        This is passed a raw TileDB object opened in read mode, since writers
        will need to retrieve data from the backing store on setup.
        """
        # non–attrs-managed field
        self.metadata = MetadataWrapper(self, dict(reader.meta))

    @property
    def reader(self) -> _RawHdl_co:
        """Accessor to assert that you are working in read mode."""
        if self.closed:
            raise SOMAError(f"{self} is closed")
        if self.mode == "r":
            return self._handle
        raise SOMAError(f"cannot read from {self}; it is open for writing")

    @property
    def writer(self) -> _RawHdl_co:
        """Accessor to assert that you are working in write mode."""
        if self.closed:
            raise SOMAError(f"{self} is closed")
        if self.mode == "w":
            return self._handle
        raise SOMAError(f"cannot write to {self}; it is open for reading")

    def close(self) -> None:
        if self.closed:
            return
        self.metadata._write()
        self._handle.close()
        self.closed = True

    def _flush_hack(self) -> None:
        """On write handles, flushes pending writes. Does nothing to reads."""
        if self.mode == "w":
            self.metadata._write()
            self._handle.close()
            self._handle = self._opener(self.uri, "w", self.context, self.timestamp_ms)

    def _check_open(self) -> None:
        if self.closed:
            raise SOMAError(f"{self!r} is closed")

    def __repr__(self) -> str:
        closed_str = " (closed)" if self.closed else ""
        return f"<{type(self).__name__} {self.mode} on {self.uri!r}{closed_str}>"

    def __enter__(self) -> Self:
        return self

    def __exit__(self, *_: Any) -> None:
        self.close()

    def __del__(self) -> None:
        self.close()


AnyWrapper = Wrapper[RawHandle]
"""Non-instantiable type representing any Handle."""


class ArrayWrapper(Wrapper[tiledb.Array]):
    """Wrapper around a TileDB Array handle."""

    @classmethod
    def _opener(
        cls,
        uri: str,
        mode: options.OpenMode,
        context: SOMATileDBContext,
        timestamp: int,
    ) -> tiledb.Array:
        return tiledb.open(
            uri,
            mode,
            timestamp=timestamp,
            ctx=context.tiledb_ctx,
        )

    @property
    def schema(self) -> tiledb.ArraySchema:
        return self._handle.schema

    def non_empty_domain(self) -> Optional[Tuple[Tuple[object, object], ...]]:
        try:
            ned = self._handle.nonempty_domain()
            if ned is None:
                return None
            return cast(Tuple[Tuple[object, object], ...], ned)
        except tiledb.TileDBError as e:
            raise SOMAError(e)

    @property
    def domain(self) -> Tuple[Tuple[object, object], ...]:
        dom = self._handle.schema.domain
        return tuple(dom.dim(i).domain for i in range(dom.ndim))

    @property
    def ndim(self) -> int:
        return int(self._handle.schema.domain.ndim)

    @property
    def attr_names(self) -> Tuple[str, ...]:
        schema = self._handle.schema
        return tuple(schema.attr(i).name for i in range(schema.nattr))

    @property
    def dim_names(self) -> Tuple[str, ...]:
        schema = self._handle.schema
        return tuple(schema.domain.dim(i).name for i in range(schema.domain.ndim))

    def enum(self, label: str) -> tiledb.Enumeration:
        return self._handle.enum(label)


@attrs.define(frozen=True)
class GroupEntry:
    uri: str
    wrapper_type: Type[AnyWrapper]

    @classmethod
    def from_object(cls, obj: tiledb.Object) -> "GroupEntry":
        if obj.type == tiledb.Array:
            return GroupEntry(obj.uri, ArrayWrapper)
        if obj.type == tiledb.Group:
            return GroupEntry(obj.uri, GroupWrapper)
        raise SOMAError(f"internal error: unknown object type {obj.type}")


class GroupWrapper(Wrapper[tiledb.Group]):
    """Wrapper around a TileDB Group handle."""

    @classmethod
    def _opener(
        cls,
        uri: str,
        mode: options.OpenMode,
        context: SOMATileDBContext,
        timestamp: int,
    ) -> tiledb.Group:
        # We want to do open-group-at-timestamp.
        ctx = context.tiledb_ctx
        cfgdict = context.tiledb_ctx.config().dict()
        cfgdict["sm.group.timestamp_end"] = timestamp
        return tiledb.Group(uri, mode, ctx=ctx, config=tiledb.Config(cfgdict))

    def _do_initial_reads(self, reader: tiledb.Group) -> None:
        super()._do_initial_reads(reader)
        self.initial_contents = {
            o.name: GroupEntry.from_object(o) for o in reader if o.name is not None
        }


class DataFrameWrapper(Wrapper[clib.SOMADataFrame]):
    """Wrapper around a Pybind11 SOMADataFrame handle."""

    @classmethod
    def _opener(
        cls,
        uri: str,
        mode: options.OpenMode,
        context: SOMATileDBContext,
        timestamp: int,
    ) -> clib.SOMADataFrame:
        open_mode = clib.OpenMode.read if mode == "r" else clib.OpenMode.write
        config = {k: str(v) for k, v in context.tiledb_config.items()}
        column_names: List[str] = []
        result_order = clib.ResultOrder.automatic
        return clib.SOMADataFrame.open(
            uri,
            open_mode,
            config,
            column_names,
            result_order,
            (0, timestamp),
        )

    # Covariant types should normally not be in parameters, but this is for
    # internal use only so it's OK.
    def _do_initial_reads(self, reader: _RawHdl_co) -> None:  # type: ignore[misc]
        """Final setup step before returning the Handle.

        This is passed a raw TileDB object opened in read mode, since writers
        will need to retrieve data from the backing store on setup.
        """
        # non–attrs-managed field
        self.metadata = MetadataWrapper(self, dict(reader.meta))

    @property
    def schema(self) -> pa.Schema:
        return self._handle.schema

    @property
    def meta(self) -> "MetadataWrapper":
        return MetadataWrapper(self, dict(self._handle.meta))

    @property
    def ndim(self) -> int:
        return len(self._handle.index_column_names)

    @property
    def count(self) -> int:
        return int(self._handle.count)

    def _cast_domain(
        self, domain: Callable[[str, DTypeLike], Tuple[object, object]]
    ) -> Tuple[Tuple[object, object], ...]:
        result = []
        for name in self._handle.index_column_names:
            dtype = self._handle.schema.field(name).type
            if pa.types.is_timestamp(dtype):
                np_dtype = np.dtype(dtype.to_pandas_dtype())
                dom = domain(name, np_dtype)
                result.append(
                    (
                        np_dtype.type(dom[0], dtype.unit),
                        np_dtype.type(dom[1], dtype.unit),
                    )
                )
            else:
                if pa.types.is_large_string(dtype) or pa.types.is_string(dtype):
                    dtype = np.dtype("U")
                elif pa.types.is_large_binary(dtype) or pa.types.is_binary(dtype):
                    dtype = np.dtype("S")
                else:
                    dtype = np.dtype(dtype.to_pandas_dtype())
                result.append(domain(name, dtype))
        return tuple(result)

    @property
    def domain(self) -> Tuple[Tuple[object, object], ...]:
        return self._cast_domain(self._handle.domain)

    def non_empty_domain(self) -> Optional[Tuple[Tuple[object, object], ...]]:
        result = self._cast_domain(self._handle.non_empty_domain)
        return result or None

    @property
    def attr_names(self) -> Tuple[str, ...]:
        return tuple(
            f.name for f in self.schema if f.name not in self._handle.index_column_names
        )

    @property
    def dim_names(self) -> Tuple[str, ...]:
        return tuple(self._handle.index_column_names)

    def enum(self, label: str) -> tiledb.Enumeration:
        # The DataFrame handle may either be ArrayWrapper or DataFrameWrapper.
        # enum is only used in the DataFrame write path and is implemented by
        # ArrayWrapper. If enum is called in the read path, it is an error.
        raise NotImplementedError


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

    owner: Wrapper[RawHandle]
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
        meta = self.owner.writer.meta
        for key, mod in self._mods.items():
            if mod in (_DictMod.ADDED, _DictMod.UPDATED):
                meta[key] = self.cache[key]
            if mod is _DictMod.DELETED:
                del meta[key]
        # Temporary hack: When we flush writes, note that the cache
        # is back in sync with disk.
        self._mods.clear()

    def __repr__(self) -> str:
        prefix = f"{type(self).__name__}({self.owner})"
        if self.owner.closed:
            return f"<{prefix}>"
        return f"<{prefix} {self.cache}>"


def _check_metadata_type(key: str, obj: object) -> None:
    """Pre-checks that a metadata entry can be stored in an array.

    These checks are reproduced from the TileDB Python metadata-setting methods,
    but are slightly more restrictive than what TileDB allows in general:
    TileDB allows (some) arrays as metadata values, but the SOMA spec does not
    allow arrays of any kind.

    We have to pre-check since we don't write metadata changes until closing.
    """
    if not isinstance(key, str):
        raise TypeError(f"metadata keys must be strings, not {type(key)}")
    if isinstance(obj, (bytes, float, int, str)):
        return
    raise TypeError(f"cannot store {type(obj)} instance as metadata")
