from typing import List, Optional, Sequence, Tuple

import pyarrow as pa
import tiledb

# This package's pybind11 code
import tiledbsoma.libtiledbsoma as clib

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
        context: Optional[SOMATileDBContext] = None,
    ):
        """
        See the ``TileDBObject`` constructor.
        """
        super().__init__(uri, context=context)

    @property
    def schema(self) -> pa.Schema:
        """
        Return data schema, in the form of an Arrow Schema.
        """
        return get_arrow_schema_from_tiledb_uri(self.uri, self._ctx)

    # tiledb.Array handle for [re]use while self is "open"
    _open_tiledb_array: Optional[tiledb.Array] = None

    def _sub_open(self) -> None:
        assert self._open_mode in ("r", "w") and self._open_tiledb_array is None
        self._open_tiledb_array = self._close_stack.enter_context(
            tiledb.open(self._uri, mode=self._open_mode, ctx=self._ctx)
        )

    @property
    def _tiledb_obj(self) -> tiledb.Array:
        "get the open tiledb.Array handle"
        assert self._open_tiledb_array is not None  # => not self.closed
        return self._open_tiledb_array

    def close(self) -> None:
        self._open_tiledb_array = None
        super().close()  # closes self._open_tiledb_array via self._close_stack

    def _tiledb_array_schema(self) -> tiledb.ArraySchema:
        """
        Returns the TileDB array schema. Not part of the SOMA API; for dev/debug/etc.
        """
        with self._ensure_open():
            return self._tiledb_obj.schema

    def _tiledb_array_keys(self) -> Sequence[str]:
        """
        Return all dim and attr names.
        """
        with self._ensure_open():
            A = self._tiledb_obj
            dim_names = [A.domain.dim(i).name for i in range(A.domain.ndim)]
            attr_names = [A.schema.attr(i).name for i in range(A.schema.nattr)]
            return dim_names + attr_names

    def _tiledb_dim_names(self) -> Tuple[str, ...]:
        """
        Reads the dimension names from the schema: for example, ['obs_id', 'var_id'].
        """
        with self._ensure_open():
            A = self._tiledb_obj
            return tuple([A.domain.dim(i).name for i in range(A.domain.ndim)])

    def _tiledb_attr_names(self) -> List[str]:
        """
        Reads the attribute names from the schema: for example, the list of column names in a dataframe.
        """
        with self._ensure_open():
            A = self._tiledb_obj
            return [A.schema.attr(i).name for i in range(A.schema.nattr)]

    def _soma_reader(
        self,
        schema: Optional[tiledb.ArraySchema] = None,
        column_names: Optional[Sequence[str]] = None,
        query_condition: Optional[tiledb.QueryCondition] = None,
        result_order: Optional[str] = None,
    ) -> clib.SOMAReader:
        """
        Construct a C++ SOMAReader using appropriate context/config/etc.
        """
        kwargs = {
            "name": self.__class__.__name__,
            "platform_config": self._ctx.config().dict(),
        }
        # Leave empty arguments out of kwargs to allow C++ constructor defaults to apply, as
        # they're not all wrapped in std::optional<>.
        if schema:
            kwargs["schema"] = schema
        if column_names:
            kwargs["column_names"] = column_names
        if query_condition:
            kwargs["query_condition"] = query_condition
        if result_order:
            kwargs["result_order"] = result_order
        return clib.SOMAReader(self._uri, **kwargs)
