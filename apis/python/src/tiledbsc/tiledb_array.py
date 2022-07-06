from typing import Dict, Optional, Sequence

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
        self, uri: str, name: str, *, parent: Optional["tiledbsc.TileDBGroup"] = None
    ):
        """
        See the TileDBObject constructor.
        """
        super().__init__(uri, name, parent=parent)

    def _open(self, mode: str = "r") -> tiledb.Array:
        """
        This is just a convenience wrapper allowing 'with self._open() as A: ...' rather than
        'with tiledb.open(self.uri) as A: ...'.
        """
        assert mode in ["w", "r"]
        # This works in either 'with self._open() as A:' or 'A = self._open(); ...; A.close().  The
        # reason is that with-as invokes our return value's __enter__ on return from this method,
        # and our return value's __exit__ on exit from the body of the with-block. The tiledb
        # array object does both of those things. (And if it didn't, we'd get a runtime AttributeError
        # on with-as, flagging the non-existence of the __enter__ or __exit__.)
        return tiledb.open(self.uri, mode=mode, ctx=self._ctx)

    def exists(self) -> bool:
        """
        Tells whether or not there is storage for the array. This might be in case a SOMA
        object has not yet been populated, e.g. before calling `from_anndata` -- or, if the
        SOMA has been populated but doesn't have this member (e.g. not all SOMAs have a `varp`).
        """
        return bool(tiledb.array_exists(self.uri))

    def tiledb_array_schema(self) -> tiledb.ArraySchema:
        """
        Returns the TileDB array schema.
        """
        with self._open() as A:
            return A.schema

    def dim_names(self) -> Sequence[str]:
        """
        Reads the dimension names from the schema: for example, ['obs_id', 'var_id'].
        """
        with self._open() as A:
            return [A.schema.domain.dim(i).name for i in range(A.schema.domain.ndim)]

    def dim_names_to_types(self) -> Dict[str, str]:
        """
        Returns a dict mapping from dimension name to dimension type.
        """
        with self._open() as A:
            dom = A.schema.domain
            return {dom.dim(i).name: dom.dim(i).dtype for i in range(dom.ndim)}

    def attr_names(self) -> Sequence[str]:
        """
        Reads the attribute names from the schema: for example, the list of column names in a dataframe.
        """
        with self._open() as A:
            return [A.schema.attr(i).name for i in range(A.schema.nattr)]

    def attr_names_to_types(self) -> Dict[str, str]:
        """
        Returns a dict mapping from attribute name to attribute type.
        """
        with self._open() as A:
            schema = A.schema
            return {
                schema.attr(i).name: schema.attr(i).dtype for i in range(schema.nattr)
            }

    def has_attr_name(self, attr_name: str) -> bool:
        """
        Returns true if the array has the specified attribute name, false otherwise.
        """
        return attr_name in self.attr_names()

    def has_attr_names(self, attr_names: Sequence[str]) -> bool:
        """
        Returns true if the array has all of the specified attribute names, false otherwise.
        """
        attr_names_set = set(self.attr_names())
        return all([attr_name in attr_names_set for attr_name in attr_names])

    def show_metadata(self, recursively: bool = True, indent: str = "") -> None:
        """
        Shows metadata for the array.
        """
        print(f"{indent}[{self.name}]")
        for key, value in self.metadata().items():
            print(f"{indent}- {key}: {value}")
