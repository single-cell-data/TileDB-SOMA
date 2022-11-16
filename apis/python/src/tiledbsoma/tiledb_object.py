from abc import ABC, abstractmethod
from typing import Any, Mapping, Optional, Sequence, Union

import tiledb

import tiledbsoma

from . import util
from .soma_options import SOMAOptions


class TileDBObject(ABC):
    """
    Base class for ``TileDBArray`` and ``TileDBGroup``.

    Manages ``soma_options``, ``ctx``, etc. which are common to both.
    """

    uri: str
    name: str
    nested_name: str
    _soma_options: SOMAOptions

    def __init__(
        self,
        # All objects:
        uri: str,
        name: str,
        *,
        # Non-top-level objects can have a parent to propgate context, depth, etc.
        parent: Optional["tiledbsoma.TileDBGroup"] = None,
        # Top-level objects should specify these:
        soma_options: Optional[SOMAOptions] = None,
        ctx: Optional[tiledb.Ctx] = None,
    ):
        """
        Initialization-handling shared between ``TileDBArray`` and ``TileDBGroup``.  Specify soma_options
        and ctx for the top-level object; omit them and specify parent for non-top-level
        objects. Note that the parent reference is solely for propagating options, ctx, display
        depth, etc.
        """
        self.uri = uri
        self.name = name

        if ctx is None:
            ctx = self._default_ctx()

        if parent is None:
            self._ctx = ctx
            self._indent = ""
            self.nested_name = name
        else:
            soma_options = parent._soma_options
            self._ctx = parent._ctx
            self._indent = parent._indent + "  "
            self.nested_name = parent.nested_name + "/" + name

        self._soma_options = soma_options or SOMAOptions()
        # Null ctx is OK if that's what they wanted (e.g. not doing any TileDB-Cloud ops).

    def _default_ctx(self) -> tiledb.Ctx:
        """
        Default TileDB configuration parameters, for when none other has been specified by the
        user.
        """
        return tiledb.Ctx(
            {
                # This is necessary for smaller tile capacities when querying with a smaller memory
                # budget.
                "sm.mem.reader.sparse_global_order.ratio_array_data": 0.3,
            }
        )

    def metadata(self) -> Mapping[str, Any]:
        """
        Returns metadata from the group/array as a dict.
        """
        with self._open("r") as obj:
            return dict(obj.meta)

    def has_metadata(self, key: str) -> bool:
        """
        Returns whether metadata is associated with the group/array.
        """
        with self._open("r") as obj:
            return key in obj.meta

    def metadata_keys(self) -> Sequence[str]:
        """
        Returns metadata keys associated with the group/array.
        """
        with self._open("r") as obj:
            return list(obj.meta.keys())

    def get_metadata(self, key: str) -> Any:
        """
        Returns metadata associated with the group/array.
        Raises ``KeyError`` if there is no such key in the metadata.
        """
        with self._open("r") as obj:
            return obj.meta[key]

    def set_metadata(self, key: str, value: Any) -> None:
        """
        Returns metadata associated with the group/array.
        """
        with self._open("w") as obj:
            obj.meta[key] = value

    def get_object_type(self) -> str:
        """
        Returns the class name associated with the group/array.
        """
        with self._open("r") as obj:
            return str(obj.meta[util.SOMA_OBJECT_TYPE_METADATA_KEY])

    def _set_object_type_metadata(self) -> None:
        """
        This helps nested-structured traversals (especially those that start at the SOMACollection
        level) confidently navigate with a minimum of introspection on group contents.
        """
        with self._open("w") as obj:
            obj.meta.update(
                {
                    util.SOMA_OBJECT_TYPE_METADATA_KEY: self.__class__.__name__,
                    util.SOMA_ENCODING_VERSION_METADATA_KEY: util.SOMA_ENCODING_VERSION,
                }
            )

    @abstractmethod
    def _open(self, mode: str = "r") -> Union[tiledb.Array, tiledb.Group]:
        """Open the underlying TileDB array or Group"""
