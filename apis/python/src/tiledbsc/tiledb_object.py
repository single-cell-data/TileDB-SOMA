from abc import ABC, abstractmethod
from typing import Any, Dict, Optional, Sequence, Union

import tiledb

import tiledbsc

from .soma_options import SOMAOptions


class TileDBObject(ABC):
    """
    Base class for `TileDBArray` and `TileDBGroup`.

    Manages soma_options, context, etc. which are common to both.
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
        parent: Optional["tiledbsc.TileDBGroup"] = None,
        # Top-level objects should specify these:
        soma_options: Optional[SOMAOptions] = None,
        ctx: Optional[tiledb.Ctx] = None,
    ):
        """
        Initialization-handling shared between `TileDBArray` and `TileDBGroup`.  Specify soma_options
        and ctx for the top-level object; omit them and specify parent for non-top-level
        objects. Note that the parent reference is solely for propagating options, ctx, display
        depth, etc.
        """
        self.uri = uri
        self.name = name

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

    def _object_type(self) -> str:
        """
        This should be implemented by child classes and should return what `tiledb.object_type(uri)`
        returns for objects of a given type -- nominally `"group"` or `"array"`.
        """
        raise Exception("This virtual method must be overridden by a child class.")

    def exists(self) -> bool:
        found = tiledb.object_type(self.uri, ctx=self._ctx)
        if found is None:
            return False
        elif found == self._object_type():
            return True
        else:
            raise Exception(
                f"Internal error: expected _object_type {self._object_type()} but found {found} at {self.uri}."
            )

    def metadata(self) -> Dict:
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
        Raises `KeyError` if there is no such key in the metadata.
        """
        with self._open("r") as obj:
            return obj.meta[key]

    def set_metadata(self, key: str, value: Any) -> None:
        """
        Returns metadata associated with the group/array.
        """
        with self._open("w") as obj:
            obj.meta[key] = value

    @abstractmethod
    def _open(self, mode: str = "r") -> Union[tiledb.Array, tiledb.Group]:
        """Open the underlying TileDB array or Group"""
