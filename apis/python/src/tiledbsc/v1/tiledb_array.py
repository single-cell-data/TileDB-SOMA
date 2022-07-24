from typing import Optional

import tiledb

import tiledbsc

from .tiledb_object import TileDBObject


class TileDBArray(TileDBObject):
    """
    Wraps arrays from TileDB-Py by retaining a URI, options, etc.
    Also serves as an abstraction layer to hide TileDB-specific details from the API, unless
    requested.
    """

    def __init__(
        self,
        uri: str,
        *,
        name: Optional[str] = None,
        parent: Optional["tiledbsc.v1.SOMACollection"] = None,
    ):
        """
        See the TileDBObject constructor.
        """
        super().__init__(uri, name=name, parent=parent)

    # TODO
    #    def delete(uri: str) -> None
    #        """
    #        Delete the SOMADataFrame specified with the URI.
    #        """

    def _tiledb_open(self, mode: str = "r") -> tiledb.Array:
        """
        This is just a convenience wrapper allowing 'with self._tiledb_open() as A: ...' rather than
        'with tiledb.open(self._uri) as A: ...'.
        """
        assert mode in ["w", "r"]
        # This works in either 'with self._tiledb_open() as A:' or 'A = self._tiledb_open(); ...; A.close().  The
        # reason is that with-as invokes our return value's __enter__ on return from this method,
        # and our return value's __exit__ on exit from the body of the with-block. The tiledb
        # array object does both of those things. (And if it didn't, we'd get a runtime AttributeError
        # on with-as, flagging the non-existence of the __enter__ or __exit__.)
        return tiledb.open(self._uri, mode=mode, ctx=self._ctx)

    def _tiledb_array_schema(self) -> tiledb.ArraySchema:
        """
        Returns the TileDB array schema. Not part of the SOMA API; for dev/debug/etc.
        """
        with self._tiledb_open() as A:
            return A.schema

    def _show_metadata(self, recursively: bool = True, indent: str = "") -> None:
        """
        Shows metadata for the array.
        XXX TEMP TEMP TEMP
        """
        print(f"{indent}[{self._name}]")
        for key, value in self.metadata:
            print(f"{indent}- {key}: {value}")
