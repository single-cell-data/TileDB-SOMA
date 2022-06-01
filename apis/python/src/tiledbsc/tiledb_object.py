import tiledb
from .soma_options import SOMAOptions

from typing import Optional

import os


class TileDBObject:
    """
    Base class for `TileDBArray` and `TileDBGroup`. Manages soma_options, context, etc. which are common
    to both.
    """

    uri: str
    name: str

    _soma_options: SOMAOptions
    _verbose: bool
    _ctx: Optional[tiledb.Ctx]

    _indent: str  # for display strings

    def __init__(
        self,
        # All objects:
        uri: str,
        name: str,
        # Non-top-level objects can have a parent to propgate context, depth, etc.
        # Circular import if we say this, but it must be a TileDBGroup:
        # parent: Optional[TileDBGroup] = None,
        parent=None,
        # Top-level objects should specify these:
        soma_options: Optional[SOMAOptions] = None,
        verbose: Optional[bool] = True,
        ctx: Optional[tiledb.Ctx] = None,
    ):
        """
        Initialization-handling shared between `TileDBArray` and `TileDBGroup`.  Specify soma_options,
        verbose, and ctx for the top-level object; omit them and specify parent for non-top-level
        objects. Note that the parent reference is solely for propagating options, ctx, display
        depth, etc.
        """

        self.uri = uri
        self.name = name

        if parent is None:
            self._soma_options = soma_options
            self._verbose = verbose
            self._ctx = ctx
            self._indent = ""
        else:
            self._soma_options = parent._soma_options
            self._verbose = parent._verbose
            self._ctx = parent._ctx
            self._indent = parent._indent + "  "

        if os.getenv("TILEDBSC_PY_SUPPRESS_VERBOSE") != None:
            self._verbose = False

        if self._soma_options is None:
            self._soma_options = SOMAOptions()
        # Null ctx is OK if that's what they wanted (e.g. not doing any TileDB-Cloud ops).

    def _object_type(self):
        """
        This should be implemented by child classes and should return what `tiledb.object_type(uri)`
        returns for objects of a given type -- nominally `"group"` or `"array"`.
        """
        raise Exception("This virtual method must be overridden by a child class.")

    def exists(self):
        found = tiledb.object_type(self.uri)
        if found == None:
            return False
        elif found == self._object_type():
            return True
        else:
            raise Exception(
                f"Internal error: expected _object_type {self._object_type()} but found {found} at {self.uri}."
            )
