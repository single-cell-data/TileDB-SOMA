import os
from abc import ABC, abstractmethod
from typing import Optional, Union

import tiledb

import tiledbsc

from . import util
from .soma_metadata_mapping import SOMAMetadataMapping
from .tiledb_platform_config import TileDBPlatformConfig


class TileDBObject(ABC):
    """
    Base class for `TileDBArray` and `SOMACollection`.

    Manages tiledb_platform_config, context, etc. which are common to both.
    """

    _uri: str
    _name: str
    _nested_name: str
    _tiledb_platform_config: TileDBPlatformConfig
    metadata: SOMAMetadataMapping

    def __init__(
        self,
        # All objects:
        uri: str,
        name: Optional[str] = None,
        *,
        # Non-top-level objects can have a parent to propgate context, depth, etc.
        parent: Optional["tiledbsc.v1.SOMACollection"] = None,
        # Top-level objects should specify these:
        tiledb_platform_config: Optional[TileDBPlatformConfig] = None,
        ctx: Optional[tiledb.Ctx] = None,
    ):
        """
        Initialization-handling shared between `TileDBArray` and `SOMACollection`.  Specify tiledb_platform_config
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

        self.metadata = SOMAMetadataMapping(self)

    def __repr__(self) -> str:
        """
        XXX TEMP
        """
        return f"name={self._name},uri={self._uri}"

    @abstractmethod
    def _tiledb_open(self, mode: str = "r") -> Union[tiledb.Array, tiledb.Group]:
        """Open the underlying TileDB array or Group"""

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

    def _get_object_type_from_metadata(self) -> str:
        """
        Returns the class name associated with the group/array.
        """
        return self.metadata.get(util.SOMA_OBJECT_TYPE_METADATA_KEY)

    def _set_object_type_metadata(self) -> None:
        """
        This helps nested-structured traversals (especially those that start at the SOMACollection
        level) confidently navigate with a minimum of introspection on group contents.
        """
        # TODO: make a multi-set in SOMAMetadataMapping that would above a double-open there.
        with self._tiledb_open("w") as obj:
            obj.meta.update(
                {
                    util.SOMA_OBJECT_TYPE_METADATA_KEY: self.__class__.__name__,
                    util.SOMA_ENCODING_VERSION_METADATA_KEY: util.SOMA_ENCODING_VERSION,
                }
            )
