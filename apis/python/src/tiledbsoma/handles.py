from typing import (
    Any,
    Callable,
    Dict,
    Generic,
    Iterator,
    MutableMapping,
    Optional,
    Set,
    TypeVar,
    Union,
)

import attrs
import tiledb
from somacore import options
from typing_extensions import Literal

from .exception import DoesNotExistError, SOMAError
from .options import SOMATileDBContext
from .util_tiledb import is_does_not_exist_error

StorageType = Literal["array", "group"]
"""How a SOMA object is stored."""
TDBHandle = Union[tiledb.Array, tiledb.Group]
"""Handle on some persistent TileDB object."""


_TDBO = TypeVar("_TDBO", bound=TDBHandle)
_TDBO_co = TypeVar("_TDBO_co", bound=TDBHandle, covariant=True)
_Self = TypeVar("_Self")


# TODO: We should only need to keep around one of the read or write handles
# after we do the initial open (where we read metadata and group contents).
# For now, we're leaving this as-is.
@attrs.define(slots=False)
class ReadWriteHandle(Generic[_TDBO_co]):
    """Paired read/write handles to a TileDB object.

    You can (and should) access the members, but do not reassign them.
    """

    uri: str
    context: SOMATileDBContext
    reader: _TDBO_co
    maybe_writer: Optional[_TDBO_co]
    closed: bool = False

    @classmethod
    def open_array(
        cls,
        uri: str,
        mode: options.OpenMode,
        context: SOMATileDBContext,
    ) -> "ReadWriteHandle[tiledb.Array]":
        return cls._open(tiledb.open, uri, mode, context)

    @classmethod
    def open_group(
        cls,
        uri: str,
        mode: options.OpenMode,
        context: SOMATileDBContext,
    ) -> "ReadWriteHandle[tiledb.Group]":
        return cls._open(tiledb.Group, uri, mode, context)

    @classmethod
    def _open(
        cls,
        tiledb_cls: Callable[..., _TDBO],
        uri: str,
        mode: options.OpenMode,
        context: SOMATileDBContext,
    ) -> "ReadWriteHandle[_TDBO]":
        try:
            rd = tiledb_cls(uri, "r", ctx=context.tiledb_ctx)
            try:
                wr = (
                    tiledb_cls(uri, "w", ctx=context.tiledb_ctx)
                    if mode == "w"
                    else None
                )
            except Exception:
                rd.close()
                raise
        except tiledb.TileDBError as tdbe:
            if is_does_not_exist_error(tdbe):
                raise DoesNotExistError(
                    f"{tiledb_cls.__name__} {uri!r} does not exist"
                ) from tdbe
            raise
        return ReadWriteHandle(uri, context, reader=rd, maybe_writer=wr)

    def __attrs_post_init__(self) -> None:
        # non–attrs-managed field
        self.metadata = MetadataWrapper(self, dict(self.reader.meta))

    @property
    def writer(self) -> _TDBO_co:
        if self.maybe_writer is None:
            raise SOMAError(f"cannot write to {self.uri!r} opened in read-only mode")
        return self.maybe_writer

    @property
    def mode(self) -> options.OpenMode:
        return "r" if self.maybe_writer is None else "w"

    def close(self) -> None:
        if self.closed:
            return
        self.closed = True
        self.reader.close()
        if self.maybe_writer is not None:
            self.maybe_writer.close()

    def _flush_hack(self) -> None:
        """Flushes writes. Use rarely, if ever."""
        if self.maybe_writer is None:
            return
        self.maybe_writer.close()
        if isinstance(self.maybe_writer, tiledb.Array):
            self.maybe_writer = tiledb.open(self.uri, "w", ctx=self.context.tiledb_ctx)
            self.reader.reopen()
        else:
            self.maybe_writer = tiledb.Group(self.uri, "w", ctx=self.context.tiledb_ctx)
            self.reader = tiledb.Group(self.uri, "r", ctx=self.context.tiledb_ctx)

    def _check_open(self) -> None:
        if self.closed:
            raise SOMAError(f"{self!r} is closed")

    def __enter__(self: _Self) -> _Self:
        return self

    def __exit__(self, *_: Any) -> None:
        self.close()

    def __del__(self) -> None:
        self.close()


@attrs.define(frozen=True)
class MetadataWrapper(MutableMapping[str, Any]):
    """A wrapper storing the metadata of some TileDB object.

    Because the view of metadata does not change after open time, we immediately
    cache all of it and use that to handle all reads. Writes are then proxied
    through to the backing store and the cache is updated to match.
    """

    owner: ReadWriteHandle[TDBHandle]
    cache: Dict[str, Any]
    _mutated_fields: Set[str] = attrs.field(factory=set)
    """Fields that have been mutated. """

    def __len__(self) -> int:
        return len(self.cache)

    def __iter__(self) -> Iterator[str]:
        return iter(self.cache)

    def __getitem__(self, key: str) -> Any:
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
