from typing import Optional, Sequence, Tuple, Type, TypeVar

import pyarrow as pa
import tiledb
from somacore import options

# This package's pybind11 code
import tiledbsoma.libtiledbsoma as clib

from .options.soma_tiledb_context import SOMATileDBContext
from .tiledb_object import TileDBObject
from .util_arrow import tiledb_schema_to_arrow
from .util_tiledb import ReadWriteHandle

_Self = TypeVar("_Self", bound="TileDBArray")


class TileDBArray(TileDBObject[tiledb.Array]):
    """
    Wraps arrays from TileDB-Py by retaining a URI, options, etc.  Also serves as an abstraction layer to hide TileDB-specific details from the API, unless requested.

    [lifecycle: experimental]
    """

    _tiledb_type = tiledb.Array

    @classmethod
    def open(
        cls: Type[_Self],
        uri: str,
        mode: options.OpenMode = "r",
        *,
        platform_config: Optional[options.PlatformConfig] = None,
        context: Optional[SOMATileDBContext] = None,
    ) -> _Self:
        del platform_config  # unused
        context = context or SOMATileDBContext()
        handle = ReadWriteHandle.open_array(uri, mode, context)

        return cls(
            handle,
            _dont_call_this_use_create_or_open_instead="tiledbsoma-internal-code",
        )

    _STORAGE_TYPE = "array"

    @property
    def schema(self) -> pa.Schema:
        """
        Return data schema, in the form of an Arrow Schema.
        """
        return tiledb_schema_to_arrow(self._tiledb_array_schema())

    def _tiledb_array_schema(self) -> tiledb.ArraySchema:
        """
        Returns the TileDB array schema, for internal use.
        """
        return self._handle.reader.schema

    def _tiledb_array_keys(self) -> Tuple[str, ...]:
        """
        Return all dim and attr names.
        """
        return self._tiledb_dim_names() + self._tiledb_attr_names()

    def _tiledb_dim_names(self) -> Tuple[str, ...]:
        """
        Reads the dimension names from the schema: for example, ['obs_id', 'var_id'].
        """
        arr = self._handle.reader
        return tuple(arr.domain.dim(i).name for i in range(arr.domain.ndim))

    def _tiledb_attr_names(self) -> Tuple[str, ...]:
        """
        Reads the attribute names from the schema: for example, the list of column names in a dataframe.
        """
        arr = self._handle.reader
        return tuple(arr.schema.attr(i).name for i in range(arr.schema.nattr))

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
        return clib.SOMAReader(self.uri, **kwargs)

    @classmethod
    def _create_internal(
        cls, uri: str, schema: tiledb.ArraySchema, context: SOMATileDBContext
    ) -> ReadWriteHandle[tiledb.Array]:
        """Creates the TileDB Array for this type and returns an opened handle.

        This does the work of creating a TileDB Array with the provided schema
        at the given URI, sets the necessary metadata, and returns a handle to
        the newly-created array, open for writing.
        """
        tiledb.Array.create(uri, schema, ctx=context.tiledb_ctx)
        handle = ReadWriteHandle.open_array(uri, "w", context)
        cls._set_create_metadata(handle)
        return handle
