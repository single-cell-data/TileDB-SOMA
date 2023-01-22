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

    # lazy tiledb.Array handle for reuse while self is "open"
    _open_tiledb_array: Optional[tiledb.Array] = None

    @property
    def _tiledb_array(self) -> tiledb.Array:
        "get the open tiledb.Array handle (opening it if needed)"
        assert self._open_mode in ("r", "w")
        if self._open_tiledb_array is None:
            self._open_tiledb_array = self._close_stack.enter_context(
                tiledb.open(self._uri, mode=self._open_mode, ctx=self._ctx)
            )
        return self._open_tiledb_array

    @property
    def _tiledb_object(self) -> tiledb.Array:
        return self._tiledb_array

    def close(self) -> None:
        self._open_tiledb_array = None
        super().close()  # closes self._open_tiledb_array via self._close_stack

    def _tiledb_array_schema(self) -> tiledb.ArraySchema:
        """
        Returns the TileDB array schema. Not part of the SOMA API; for dev/debug/etc.
        """
        with self._maybe_open():
            return self._tiledb_array.schema

    def _tiledb_array_keys(self) -> Sequence[str]:
        """
        Return all dim and attr names.
        """
        with self._maybe_open():
            A = self._tiledb_array
            dim_names = [A.domain.dim(i).name for i in range(A.domain.ndim)]
            attr_names = [A.schema.attr(i).name for i in range(A.schema.nattr)]
            return dim_names + attr_names

    def _tiledb_dim_names(self) -> Tuple[str, ...]:
        """
        Reads the dimension names from the schema: for example, ['obs_id', 'var_id'].
        """
        with self._maybe_open():
            A = self._tiledb_array
            return tuple([A.domain.dim(i).name for i in range(A.domain.ndim)])

    def _tiledb_attr_names(self) -> List[str]:
        """
        Reads the attribute names from the schema: for example, the list of column names in a dataframe.
        """
        with self._maybe_open():
            A = self._tiledb_array
            return [A.schema.attr(i).name for i in range(A.schema.nattr)]
