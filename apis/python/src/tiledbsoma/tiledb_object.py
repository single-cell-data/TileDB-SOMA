from abc import ABC, abstractproperty
from contextlib import ExitStack, contextmanager
from types import TracebackType
from typing import Iterator, Optional, Type, TypeVar, Union

import somacore
import tiledb

from . import util
from .metadata_mapping import MetadataMapping
from .options import SOMATileDBContext

# type variable for methods returning self; see PEP 673
TTileDBObject = TypeVar("TTileDBObject", bound="TileDBObject")


# TODO: SOMAObject should inherit typing.ContextManager
class TileDBObject(ABC, somacore.SOMAObject):
    """
    Base class for ``TileDBArray`` and ``Collection``.

    Accepts a SOMATileDBContext, to enable session state to be shared across SOMA objects.
    """

    _uri: str
    _context: SOMATileDBContext
    _metadata: MetadataMapping

    # "r", "w", or falsey
    _open_mode: str = ""
    # iff open: child contexts to exit when self is exited/closed
    _close_stack: ExitStack

    def __init__(
        self,
        # All objects:
        uri: str,
        *,
        # Non-top-level objects can have a parent to propagate context, depth, etc.
        parent: Optional["TileDBObject"] = None,
        # Top-level objects should specify this:
        context: Optional[SOMATileDBContext] = None,
    ):
        """
        Initialization-handling shared between ``TileDBArray`` and ``Collection``.  Specify ``context`` for
        the top-level object; omit it and specify parent for non-top-level objects. Note that the parent reference
        is solely for propagating the context
        """

        self._uri = uri

        if parent is not None:
            if context is not None:
                raise TypeError(
                    "Only one of `context` and `parent` params can be passed as an arg"
                )
            # inherit from parent
            self._context = parent._context
        else:
            self._context = context or SOMATileDBContext()

        self._metadata = MetadataMapping(self)
        self._close_stack = ExitStack()

    @property
    def context(self) -> SOMATileDBContext:
        return self._context

    @property
    def _ctx(self) -> tiledb.Ctx:
        return self._context.tiledb_ctx

    @property
    def metadata(self) -> MetadataMapping:
        """Metadata accessor"""
        return self._metadata
        # Note: this seems trivial, like we could just have `metadata` as an attribute.
        # However, we've found that since in `somacore` it's implemented as `@property`,
        # to avoid a static-analysis failure we have to do the same here.

    def delete(self) -> None:
        """
        Delete the storage specified with the URI.

        TODO: should this raise an error if the object does not exist?
        """
        try:
            tiledb.remove(self._uri)
        except tiledb.TileDBError:
            pass
        return

    def __repr__(self) -> str:
        """
        Default repr
        """
        if self.exists():
            return f'{self.soma_type}(uri="{self._uri}")'
        else:
            return f"{self.soma_type}(not created)"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, TileDBObject):
            return False
        return self._uri == other._uri

    @property
    def uri(self) -> str:
        """
        Accessor for the object's storage URI
        """
        return self._uri

    def exists(self) -> bool:
        """
        Returns true if the object exists and has the desired class name.

        This might be in case an object has not yet been populated, or, if a containing object has been populated but doesn't have a particular member (e.g. not all ``Measurement`` objects have a ``varp``).

        For ``tiledb://`` URIs this is a REST-server request which we'd like to cache.  However, remove-and-replace use-cases are possible and common in notebooks and it turns out caching the existence-check isn't a robust approach.
        """

        # Pre-checking if the group exists by calling tiledb.object_type is simple, however, for
        # tiledb-cloud URIs that occurs a penalty of two HTTP requests to the REST server, even
        # before a third, successful HTTP request for group-open.  Instead, we directly attempt the
        # group-open request, checking for an exception.
        try:
            with self._maybe_open():
                return (
                    self._metadata.get(util.SOMA_OBJECT_TYPE_METADATA_KEY)
                    == self.soma_type
                )
        except tiledb.cc.TileDBError:
            return False

    def open(self: TTileDBObject, mode: str = "r") -> TTileDBObject:
        """
        Open the object for a series of read or write operations. This is optional, but makes a
        series of small operations more efficient by allowing reuse of handles and metadata. If
        opened, the object should later be closed to release such resources. This can be automated
        with a context manager, e.g.:

          with tiledbsoma.SparseNDArray(uri=...).open() as X:
              data = X.read_table(...)
              data = X.read_table(...)

        If the object is not opened, then any individual read or write operation effectively opens
        and closes it transiently for the operation.
        """
        if mode not in ("r", "w"):
            raise ValueError(
                self.__class__.__name__ + " open mode must be one of 'r', 'w'"
            )
        if self._open_mode:
            raise RuntimeError(self.__class__.__name__ + " is already open")
        self._open_mode = mode
        return self

    def close(self) -> None:
        """
        Release any resources held while the object is open. Closing an already-closed object is a
        no-op.
        """
        if self._open_mode:
            self._open_mode = ""
            self._close_stack.close()

    def __enter__(self: TTileDBObject) -> TTileDBObject:
        assert (
            self._open_mode
        ), f"use {self.__class__.__name__}.open() as context manager"
        return self

    def __exit__(
        self,
        exc_type: Optional[Type[BaseException]],
        exc_value: Optional[BaseException],
        traceback: Optional[TracebackType],
    ) -> None:
        if self._open_mode:
            self.close()

    @contextmanager
    def _maybe_open(self, mode: str = "r") -> Iterator[None]:
        """
        Internal helper context for read/write methods to open self just for the duration of the
        operation -- IFF self isn't already open. Also rejects attempts to write when self is
        already open read-only.

        when mode == "w":
            if self is already open in read-only mode, raise an error.
            if self is already open in w mode, do nothing.
            otherwise open("w") until context exit.
        when mode == "r":
            if self is already open (IN EITHER MODE*), do nothing.
            otherwise open("r") until context exit.

            * We don't necessarily reject reads while we're open for writing, and different
              timestamps may be used in this case. Usually, reads while open for writing just
              pertain to metadata (e.g. exists()). The underlying TileDB object may reject other
              cases.
        """
        assert mode in ("r", "w")
        if self._open_mode:
            if self._open_mode == "r" and mode == "w":
                raise RuntimeError(
                    self.__class__.__name__ + " is already open read-only"
                )
            yield
        else:
            with self.open(mode):
                yield

    def __del__(self) -> None:
        self.close()

    @abstractproperty
    def _tiledb_obj(self) -> Union[tiledb.Array, tiledb.Group]:
        """
        Get reference to open TileDB object handle (self must be open)
        """
        ...

    def _common_create(self, soma_type: str) -> None:
        """
        This helps nested-structure traversals (especially those that start at the
        Collection level) confidently navigate with a minimum of introspection on group
        contents.
        """
        with self._maybe_open("w"):
            self._tiledb_obj.meta.update(
                {
                    util.SOMA_OBJECT_TYPE_METADATA_KEY: soma_type,
                    util.SOMA_ENCODING_VERSION_METADATA_KEY: util.SOMA_ENCODING_VERSION,
                }
            )
