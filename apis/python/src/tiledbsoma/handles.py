from typing import (
    Any,
    Callable,
    Dict,
    Generic,
    Iterator,
    MutableMapping,
    Optional,
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
@attrs.define()
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

    @property
    def writer(self) -> _TDBO_co:
        if self.maybe_writer is None:
            raise SOMAError(f"cannot write to {self.uri!r} opened in read-only mode")
        return self.maybe_writer

    @property
    def mode(self) -> options.OpenMode:
        return "r" if self.maybe_writer is None else "w"

    @property
    def metadata(self) -> MutableMapping[str, Any]:
        if self.maybe_writer is not None:
            return self.maybe_writer.meta  # type: ignore[no-any-return]
        return self.reader.meta  # type: ignore[no-any-return]

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

    def __enter__(self: _Self) -> _Self:
        return self

    def __exit__(self, *_: Any) -> None:
        self.close()

    def __del__(self) -> None:
        self.close()


class MetadataMapping(MutableMapping[str, Any]):
    def __init__(self, underlying: ReadWriteHandle[tiledb.Object]):
        self._underlying = underlying

    def __delitem__(self, key: str) -> None:
        """
        Remove the key from the metadata collection.
        """
        del self._underlying.writer.meta[key]

    def __iter__(self) -> Iterator[str]:
        """
        Return an iterator over the metadata collection.
        """
        return iter(self.as_dict())

    def __len__(self) -> int:
        """
        Return the number of elements in the metadata collection.

        Returns
        -------
        int
            The number of elements in the metadata collection.
        """
        return len(self.as_dict())

    def as_dict(self) -> Dict[str, Any]:
        """
        Return the metadata collection as a ``dict``.

        Returns
        -------
        dict[str, any]
            The contents of the metadata collection.
        """
        return dict(self._underlying.reader.meta)

    def __getitem__(self, key: str) -> Any:
        """
        Return the metadata element specified by ``key``.

        Parameters
        ----------
        key : str
            The name of the item.
        """
        return self._underlying.reader.meta[key]

    def __setitem__(self, key: str, value: Any) -> None:
        """
        Update the metadata collection with a new element.

        Parameters
        ----------
        key : str
            The metadata element name.
        value : any
            The metadata element value. Must be a primitive type (int,
            float, bool) or string.

        Returns
        -------
        None
        """

        if type(key) != str:
            raise TypeError("Metadata keys must be a string.")

        self._underlying.writer.meta[key] = value
