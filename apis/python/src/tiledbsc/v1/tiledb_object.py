import os
from abc import ABC, abstractmethod
from typing import Any, Mapping, Optional, Sequence, Union

import tiledb

import tiledbsc

from . import util
from .tiledb_platform_config import TileDBPlatformConfig


class TileDBObject(ABC):
    """
    Base class for `TileDBArray` and `TileDBGroup`.

    Manages tiledb_platform_config, context, etc. which are common to both.
    """

    _uri: str
    _name: str
    _nested_name: str
    _tiledb_platform_config: TileDBPlatformConfig

    def __init__(
        self,
        # All objects:
        uri: str,
        name: Optional[str] = None,
        *,
        # Non-top-level objects can have a parent to propgate context, depth, etc.
        parent: Optional["tiledbsc.v1.TileDBGroup"] = None,
        # Top-level objects should specify these:
        tiledb_platform_config: Optional[TileDBPlatformConfig] = None,
        ctx: Optional[tiledb.Ctx] = None,
    ):
        """
        Initialization-handling shared between `TileDBArray` and `TileDBGroup`.  Specify tiledb_platform_config
        and ctx for the top-level object; omit them and specify parent for non-top-level
        objects. Note that the parent reference is solely for propagating options, ctx, display
        depth, etc.
        """
        self._uri = uri
        if name is None:
            self._name = os.path.basename(uri)
        else:
            self._name = name

        if parent is None:
            self._ctx = ctx
            self._indent = ""
            self._nested_name = self._name
        else:
            tiledb_platform_config = parent._tiledb_platform_config
            self._ctx = parent._ctx
            self._indent = parent._indent + "  "
            self._nested_name = parent._nested_name + "/" + self._name

        self._tiledb_platform_config = tiledb_platform_config or TileDBPlatformConfig()
        # Null ctx is OK if that's what they wanted (e.g. not doing any TileDB-Cloud ops).

    def __repr__(self):
        """
        XXX TEMP
        """
        return f"name={self._name},uri={self._uri}"

    @abstractmethod
    def _tiledb_open(self, mode: str = "r") -> Union[tiledb.Array, tiledb.Group]:
        """Open the underlying TileDB array or Group"""

    #    def metadata(self) -> Mapping[str, Any]:
    #        """
    #        Returns metadata from the group/array as a dict.
    #        """
    #        with self._tiledb_open("r") as obj:
    #            return dict(obj.meta)

    #    def has_metadata(self, key: str) -> bool:
    #        """
    #        Returns whether metadata is associated with the group/array.
    #        """
    #        with self._tiledb_open("r") as obj:
    #            return key in obj.meta

    #    def metadata_keys(self) -> Sequence[str]:
    #        """
    #        Returns metadata keys associated with the group/array.
    #        """
    #        with self._tiledb_open("r") as obj:
    #            return list(obj.meta.keys())

    #    def get_metadata(self, key: str) -> Any:
    #        """
    #        Returns metadata associated with the group/array.
    #        Raises `KeyError` if there is no such key in the metadata.
    #        """
    #        with self._tiledb_open("r") as obj:
    #            return obj.meta[key]

    #    def set_metadata(self, key: str, value: Any) -> None:
    #        """
    #        Returns metadata associated with the group/array.
    #        """
    #        with self._tiledb_open("w") as obj:
    #            obj.meta[key] = value

    def get_name(self) -> str:
        return self._name

    def get_uri(self) -> str:
        return self._uri

    def get_type(self) -> str:
        return type(self).__name__

    #    def get_object_type(self) -> str:
    #        """
    #        Returns the class name associated with the group/array.
    #        """
    #        with self._tiledb_open("r") as obj:
    #            return str(obj.meta[util.SOMA_OBJECT_TYPE_METADATA_KEY])

    #    def _set_object_type_metadata(self) -> None:
    #        """
    #        This helps nested-structured traversals (especially those that start at the SOMACollection
    #        level) confidently navigate with a minimum of introspection on group contents.
    #        """
    #        with self._tiledb_open("w") as obj:
    #            obj.meta.update(
    #                {
    #                    util.SOMA_OBJECT_TYPE_METADATA_KEY: self.__class__.__name__,
    #                    util.SOMA_ENCODING_VERSION_METADATA_KEY: util.SOMA_ENCODING_VERSION,
    #                }
    #            )

    # ================================================================
    # ================================================================
    # ================================================================
    # XXX TEMP WIP COPY/FACTOR FROM V0

    def _common_create(self) -> None:
        """
        TODO: COMMENT
        """
        self._set_object_type_metadata()

    def _metadata(self) -> Mapping[str, Any]:
        """
        Returns metadata from the group/array as a dict.
        """
        with self._tiledb_open("r") as obj:
            return dict(obj.meta)

    def _has_metadata(self, key: str) -> bool:
        """
        Returns whether metadata is associated with the group/array.
        """
        with self._tiledb_open("r") as obj:
            return key in obj.meta

    def _metadata_keys(self) -> Sequence[str]:
        """
        Returns metadata keys associated with the group/array.
        """
        with self._tiledb_open("r") as obj:
            return list(obj.meta.keys())

    def _get_metadata(self, key: str) -> Any:
        """
        Returns metadata associated with the group/array.
        Raises `KeyError` if there is no such key in the metadata.
        """
        with self._tiledb_open("r") as obj:
            return obj.meta[key]

    def _set_metadata(self, key: str, value: Any) -> None:
        """
        Returns metadata associated with the group/array.
        """
        with self._tiledb_open("w") as obj:
            obj.meta[key] = value

    def _get_object_type_from_metadata(self) -> str:
        """
        Returns the class name associated with the group/array.
        """
        with self._tiledb_open("r") as obj:
            return str(obj.meta[util.SOMA_OBJECT_TYPE_METADATA_KEY])

    def _set_object_type_metadata(self) -> None:
        """
        This helps nested-structured traversals (especially those that start at the SOMACollection
        level) confidently navigate with a minimum of introspection on group contents.
        """
        with self._tiledb_open("w") as obj:
            obj.meta.update(
                {
                    util.SOMA_OBJECT_TYPE_METADATA_KEY: self.__class__.__name__,
                    util.SOMA_ENCODING_VERSION_METADATA_KEY: util.SOMA_ENCODING_VERSION,
                }
            )
