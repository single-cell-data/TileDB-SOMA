from typing import List, Optional, Sequence, Tuple

import pyarrow as pa
import tiledb

from .options import SOMATileDBContext
from .tiledb_object import TileDBObject
from .util_arrow import get_arrow_schema_from_tiledb_uri


class TileDBArray(TileDBObject):
    """
    Wraps arrays from TileDB-Py by retaining a URI, options, etc.  Also serves as an abstraction layer to hide TileDB-specific details from the API, unless requested.
    """

    def __init__(
        self,
        uri: str,
        *,
        parent: Optional["TileDBObject"] = None,
        # Top-level objects should specify this:
        context: Optional[SOMATileDBContext] = None,
    ):
        """
        See the ``TileDBObject`` constructor.
        """
        super().__init__(uri, parent=parent, context=context)

    @property
    def schema(self) -> pa.Schema:
        """
        Return data schema, in the form of an Arrow Schema.
        """
        return get_arrow_schema_from_tiledb_uri(self.uri, self._ctx)

    def _tiledb_open(self, mode: str = "r") -> tiledb.Array:
        """
        This is just a convenience wrapper allowing 'with self._tiledb_open() as A: ...' rather than 'with tiledb.open(self._uri) as A: ...'.
        """
        if mode not in ["r", "w"]:
            raise ValueError(f'expected mode to be one of "r" or "w"; got "{mode}"')
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

    def _tiledb_array_keys(self) -> Sequence[str]:
        """
        Return all dim and attr names.
        """
        with self._tiledb_open() as A:
            dim_names = [A.domain.dim(i).name for i in range(A.domain.ndim)]
            attr_names = [A.schema.attr(i).name for i in range(A.schema.nattr)]
            return dim_names + attr_names

    def _tiledb_dim_names(self) -> Tuple[str, ...]:
        """
        Reads the dimension names from the schema: for example, ['obs_id', 'var_id'].
        """
        with self._tiledb_open() as A:
            return tuple([A.domain.dim(i).name for i in range(A.domain.ndim)])

    def _tiledb_attr_names(self) -> List[str]:
        """
        Reads the attribute names from the schema: for example, the list of column names in a dataframe.
        """
        with self._tiledb_open() as A:
            return [A.schema.attr(i).name for i in range(A.schema.nattr)]

    def _show_metadata(self, recursively: bool = True, indent: str = "") -> None:
        """
        Prints metadata for the array, for interactive use.
        """
        for key, value in self.metadata.items():
            print(f"{indent}- {key}: {value}")
