"""Abstractions to more easily manage read and write access to TileDB data.

``open``, ``ArrayWrapper.open``, ``GroupWrapper.open`` are the important parts.
"""

import abc
from typing import (
    Any,
    Dict,
    Generic,
    Iterator,
    MutableMapping,
    Set,
    Type,
    TypeVar,
    Union,
)

import attrs
import tiledb
from somacore import options
from typing_extensions import Self

from .exception import DoesNotExistError, SOMAError, is_does_not_exist_error
from .options import SOMATileDBContext

RawHandle = Union[tiledb.Array, tiledb.Group]
_RawHdl_co = TypeVar("_RawHdl_co", bound=RawHandle, covariant=True)
"""A raw TileDB object. Covariant because Handles are immutable enough."""


def open(
    uri: str, mode: options.OpenMode, context: SOMATileDBContext
) -> "Wrapper[RawHandle]":
    """Determine whether the URI is an array or group, and open it."""
    obj_type = tiledb.object_type(uri, ctx=context.tiledb_ctx)
    if not obj_type:
        raise DoesNotExistError(f"{uri!r} does not exist")
    if obj_type == "array":
        return ArrayWrapper.open(uri, mode, context)
    if obj_type == "group":
        return GroupWrapper.open(uri, mode, context)
    raise SOMAError(f"{uri!r} has unknown storage type {obj_type!r}")


@attrs.define(eq=False, hash=False, slots=False)
class Wrapper(Generic[_RawHdl_co], metaclass=abc.ABCMeta):
    """Wrapper for TileDB handles to manage lifecycle and metadata.

    Callers may read and use (non-underscored) members but should never set
    attributes on instances.
    """

    uri: str
    mode: options.OpenMode
    context: SOMATileDBContext
    _handle: _RawHdl_co
    closed: bool = attrs.field(default=False, init=False)

    @classmethod
    def open(cls, uri: str, mode: options.OpenMode, context: SOMATileDBContext) -> Self:
        if mode not in ("r", "w"):
            raise ValueError(f"Invalid open mode {mode!r}")
        try:
            tdb = cls._opener(uri, mode, context)
            handle = cls(uri, mode, context, tdb)
            if mode == "w":
                # Briefly open a read-mode handle to populate metadata/schema/group contents/etc.,
                # ignoring any read_timestamp set in the context to get an up-to-date view in
                # preparation for writing.
                with cls._opener(
                    uri, "r", context, use_latest_read_timestamp=True
                ) as auxiliary_reader:
                    handle._do_initial_reads(auxiliary_reader)
            else:
                handle._do_initial_reads(tdb)
        except tiledb.TileDBError as tdbe:
            if is_does_not_exist_error(tdbe):
                raise DoesNotExistError(f"{uri!r} does not exist") from tdbe
            raise
        return handle

    @classmethod
    @abc.abstractmethod
    def _opener(
        cls,
        uri: str,
        mode: options.OpenMode,
        context: SOMATileDBContext,
        use_latest_read_timestamp: bool = False,
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
        self._handle.close()
        self.closed = True

    def _flush_hack(self) -> None:
        """On write handles, flushes pending writes. Does nothing to reads."""
        if self.mode == "w":
            self._handle.close()
            self._handle = self._opener(self.uri, "w", self.context)

    def _check_open(self) -> None:
        if self.closed:
            raise SOMAError(f"{self!r} is closed")

    def __repr__(self) -> str:
        return f"<{type(self).__name__} {self.mode} on {self.uri!r}>"

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
        use_latest_read_timestamp: bool = False,
    ) -> tiledb.Array:
        if not use_latest_read_timestamp:
            timestamp_arg = (
                context.write_timestamp
                if mode == "w"
                else (context.read_timestamp_start, context.read_timestamp)
            )
        else:
            # array opened in write mode should initialize with latest metadata
            assert mode == "r"
            timestamp_arg = None
        return tiledb.open(uri, mode, timestamp=timestamp_arg, ctx=context.tiledb_ctx)

    @property
    def schema(self) -> tiledb.ArraySchema:
        return self._handle.schema


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
        use_latest_read_timestamp: bool = False,
    ) -> tiledb.Group:
        if not use_latest_read_timestamp:
            # As of Feb 2023, tiledb.Group() has no timestamp arg; instead its timestamps must be
            # set in the tiledb.Ctx config. SOMATileDBContext prepares the suitable tiledb.Ctx.
            ctx_arg = (
                context._group_write_tiledb_ctx
                if mode == "w"
                else context._group_read_tiledb_ctx
            )
        else:
            # Group opened in write mode should initialize with latest contents & metadata
            assert mode == "r"
            ctx_arg = context.tiledb_ctx
        return tiledb.Group(uri, mode, ctx=ctx_arg)

    def _do_initial_reads(self, reader: tiledb.Group) -> None:
        super()._do_initial_reads(reader)
        self.initial_contents = {
            o.name: GroupEntry.from_object(o) for o in reader if o.name is not None
        }


@attrs.define(frozen=True)
class MetadataWrapper(MutableMapping[str, Any]):
    """A wrapper storing the metadata of some TileDB object.

    Because the view of metadata does not change after open time, we immediately
    cache all of it and use that to handle all reads. Writes are then proxied
    through to the backing store and the cache is updated to match.
    """

    owner: Wrapper[RawHandle]
    cache: Dict[str, Any]
    _mutated_fields: Set[str] = attrs.field(factory=set)
    """Fields that have been mutated. """

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
        self.owner._check_open()
        self.owner.writer.meta[key] = value
        self.cache[key] = value
        self._mutated_fields.add(key)

    def __delitem__(self, key: str) -> None:
        # HACK: Currently, deleting a just-added metadata entry DOES NOT WORK.
        # We track mutated fields to ensure that this works as expected,
        # at the expense of otherwise-unnecessary flushes.
        self.owner._check_open()
        if key in self._mutated_fields:
            self.owner._flush_hack()
        del self.owner.writer.meta[key]
        del self.cache[key]
        self._mutated_fields.add(key)

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.owner})"
