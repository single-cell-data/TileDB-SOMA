import os
import time
from abc import ABC, abstractmethod, abstractproperty
from types import TracebackType
from typing import ContextManager, Optional, Type, TypeVar, Union, cast

import tiledb

from . import util
from .metadata_mapping import MetadataMapping
from .tiledb_platform_config import TileDBPlatformConfig

# type variable for methods returning self; see PEP 673
TTileDBObject = TypeVar("TTileDBObject", bound="TileDBObject")


class TileDBObject(ContextManager["TileDBObject"], ABC):
    """
    Base class for ``TileDBArray`` and ``Collection``.

    Manages tiledb_platform_config, context, etc. which are common to both.
    """

    _uri: str
    _tiledb_platform_config: TileDBPlatformConfig
    _metadata: MetadataMapping
    _parent_timestamp: Optional[int] = None

    def __init__(
        self,
        # All objects:
        uri: str,
        *,
        # Non-top-level objects can have a parent to propgate context, depth, etc.
        parent: Optional["TileDBObject"] = None,
        # Top-level objects should specify these:
        tiledb_platform_config: Optional[TileDBPlatformConfig] = None,
        ctx: Optional[tiledb.Ctx] = None,
    ):
        """
        Initialization-handling shared between ``TileDBArray`` and ``Collection``.  Specify
        ``tiledb_platform_config`` and ``ctx`` for the top-level object; omit them and specify
        parent for non-top-level objects. Note that the parent reference is solely for propagating
        options, ctx, display depth, etc.
        """

        if ctx is None:
            ctx = self._default_ctx()
        self._uri = uri
        if parent is None:
            self._ctx = ctx
        else:
            tiledb_platform_config = parent._tiledb_platform_config
            self._ctx = parent._ctx
            self._parent_timestamp = parent._effective_timestamp()

        self._tiledb_platform_config = tiledb_platform_config or TileDBPlatformConfig()
        # Null ctx is OK if that's what they wanted (e.g. not doing any TileDB-Cloud ops).

        self._metadata = MetadataMapping(self)

    def _default_ctx(self) -> tiledb.Ctx:
        """
        The TileDB context used when no other is supplied. Must have good defaults for positive
        out-of-the-box UX.
        """

        cfg = {}

        # This is necessary for smaller tile capacities when querying with a smaller memory budget.
        cfg["sm.mem.reader.sparse_global_order.ratio_array_data"] = 0.3

        # Temp workaround pending https://app.shortcut.com/tiledb-inc/story/23827
        region = os.getenv("AWS_DEFAULT_REGION")
        if region is not None:
            cfg["vfs.s3.region"] = cast(str, region)  # type: ignore

        return tiledb.Ctx(cfg)

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

    @abstractproperty
    def soma_type(self) -> str:
        """
        Returns the SOMA object type, e.g. "SOMADataFrame".
        """
        ...

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
            return self._get_soma_type_from_metadata() == self.soma_type
        except tiledb.cc.TileDBError:
            return False

    _open_mode: Optional[str] = None
    _open_timestamp: Optional[int] = None

    def open(
        self: TTileDBObject, mode: str = "r", timestamp: Optional[int] = None
    ) -> TTileDBObject:
        """
        Open the object for timestamped read or write operations. The object should be closed when
        done. It's recommended to use the object as a context manager to automate this, e.g.:

          with tiledbsoma.SparseNDArray(uri=...).open() as X:
              data = X.read_table(...)
              data = X.read_table(...)

        ``mode='r'`` (default): opening at a given timestamp (default: now) provides snapshot
        consistency over a series of read operations on the current object and (for collections)
        accessed elements.

        ``mode='w'``: opening an object for writing provides that multiple write operations on the
        object all share the same timestamp. (To write a collection with timestamps shared amongst
        all elements, explicitly open each element with the desired timestamp.)

        Opening an object for a series of read or write operations can also make them more
        efficient, by reusing resources and metadata between operations.
        """
        if self.is_open():
            raise RuntimeError("TileDBObject is already open")
        if mode not in ("r", "w"):
            raise ValueError("unknown mode")
        self._open_mode = mode
        if timestamp is not None:
            self._open_timestamp = timestamp
        elif mode == "r" and self._parent_timestamp is not None:
            # Auto-inherit the parent timestamp only for reads, for now.
            # TODO: we'd want to inherit the parent timestamp for writes only if the parent itself
            # is currently open for writing. We don't currently store the reference to parent
            # beyond the initializer. What if parent is currently open for reading, not writing?
            # What if parent isn't open but grandparent is? All these complications while the
            # typical realistic use case is that the array is created and written first, and only
            # thereafter added to a parent collection. Punting for now...
            self._open_timestamp = self._parent_timestamp
        else:
            self._open_timestamp = int(time.time() * 1000)
        return self

    def is_open(self) -> bool:
        return self._open_mode is not None

    def _effective_timestamp(self) -> Optional[int]:
        if self.is_open():
            assert self._open_timestamp is not None
            return self._open_timestamp
        return self._parent_timestamp

    def close(self) -> None:
        """
        Release any resources held while the object is open.
        """
        if not self.is_open():
            raise RuntimeError("TileDBObject is not open")
        self._open_mode = None

    def __enter__(self: TTileDBObject) -> TTileDBObject:
        return self

    def __exit__(
        self,
        exc_type: Optional[Type[BaseException]],
        exc_value: Optional[BaseException],
        traceback: Optional[TracebackType],
    ) -> None:
        if self.is_open():
            self.close()

    @abstractmethod
    def _tiledb_open(self, mode: str = "r") -> Union[tiledb.Array, tiledb.Group]:
        """Open the underlying TileDB array or Group"""
        ...

    def _common_create(self, soma_type: str) -> None:
        """
        Utility method for various constructors.
        """
        self._set_object_type_metadata(soma_type)

    def _set_object_type_metadata(self, soma_type: str) -> None:
        """
        This helps nested-structure traversals (especially those that start at the Collection level) confidently navigate with a minimum of introspection on group contents.
        """
        # TODO: make a multi-set in MetadataMapping that would avoid a double-open there.
        with self._tiledb_open("w") as obj:
            obj.meta.update(
                {
                    util.SOMA_OBJECT_TYPE_METADATA_KEY: soma_type,
                    util.SOMA_ENCODING_VERSION_METADATA_KEY: util.SOMA_ENCODING_VERSION,
                }
            )

    def _get_soma_type_from_metadata(self) -> str:
        """
        Returns the class name associated with the group/array.
        """
        return cast(str, self._metadata.get(util.SOMA_OBJECT_TYPE_METADATA_KEY))
