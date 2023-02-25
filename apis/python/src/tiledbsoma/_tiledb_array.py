from typing import Optional, Sequence, Tuple

import pyarrow as pa
import tiledb

from . import _tdb_handles, _util

# This package's pybind11 code
from . import libtiledbsoma as clib
from ._arrow_types import tiledb_schema_to_arrow
from ._tiledb_object import TileDBObject
from ._types import is_nonstringy_sequence
from .options._soma_tiledb_context import SOMATileDBContext


class TileDBArray(TileDBObject[_tdb_handles.ArrayWrapper]):
    """
    Wraps arrays from TileDB-Py by retaining a URI, options, etc.  Also serves as an abstraction layer to hide TileDB-specific details from the API, unless requested.

    [lifecycle: experimental]
    """

    __slots__ = ()

    _wrapper_type = _tdb_handles.ArrayWrapper

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
        return self._handle.schema

    def _tiledb_array_keys(self) -> Tuple[str, ...]:
        """
        Return all dim and attr names.
        """
        return self._tiledb_dim_names() + self._tiledb_attr_names()

    def _tiledb_dim_names(self) -> Tuple[str, ...]:
        """
        Reads the dimension names from the schema: for example, ['obs_id', 'var_id'].
        """
        schema = self._handle.schema
        return tuple(schema.domain.dim(i).name for i in range(schema.domain.ndim))

    def _tiledb_attr_names(self) -> Tuple[str, ...]:
        """
        Reads the attribute names from the schema: for example, the list of column names in a dataframe.
        """
        schema = self._handle.schema
        return tuple(schema.attr(i).name for i in range(schema.nattr))

    def _soma_reader(
        self,
        *,
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
            "timestamp": self.context._timestamp_arg(),
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

    def _set_reader_coords(self, sr: clib.SOMAReader, coords: Sequence[object]) -> None:
        """Parses the given coords and sets them on the SOMA Reader."""
        if not is_nonstringy_sequence(coords):
            raise TypeError(
                f"coords type {type(coords)} must be a regular sequence,"
                " not str or bytes"
            )
        schema = self._handle.schema
        if len(coords) > schema.domain.ndim:
            raise ValueError(
                f"coords ({len(coords)} elements) must be shorter than ndim"
                f" ({schema.domain.ndim})"
            )
        for i, coord in enumerate(coords):
            dim = self._handle.schema.domain.dim(i)
            if not self._set_reader_coord(sr, i, dim, coord):
                raise TypeError(
                    f"coord type {type(coord)} for dimension {dim.name}"
                    f" (slot {i}) unsupported"
                )

    def _set_reader_coord(
        self, sr: clib.SOMAReader, dim_idx: int, dim: tiledb.Dim, coord: object
    ) -> bool:
        """Parses a single coordinate entry.

        The base implementation parses the most fundamental types shared by all
        TileDB Array types; subclasses can implement their own readers that
        handle types not recognized here.

        :return: True if successful, False if unrecognized.
        """
        del dim_idx  # Unused.
        if coord is None:
            return True  # No constraint; select all in this dimension

        if isinstance(coord, int):
            sr.set_dim_points_int64(dim.name, [coord])
            return True
        if isinstance(coord, slice):
            _util.validate_slice(coord)
            try:
                lo_hi = _util.slice_to_numeric_range(coord, dim.domain)
            except _util.NonNumericDimensionError:
                return False  # We only handle numeric dimensions here.
            if lo_hi:
                sr.set_dim_ranges_int64(dim.name, [lo_hi])
            # If `None`, coord was `slice(None)` and there is no constraint.
            return True
        return False

    @classmethod
    def _create_internal(
        cls, uri: str, schema: tiledb.ArraySchema, context: SOMATileDBContext
    ) -> _tdb_handles.ArrayWrapper:
        """Creates the TileDB Array for this type and returns an opened handle.

        This does the work of creating a TileDB Array with the provided schema
        at the given URI, sets the necessary metadata, and returns a handle to
        the newly-created array, open for writing.
        """
        tiledb.Array.create(uri, schema, ctx=context.tiledb_ctx)
        handle = cls._wrapper_type.open(uri, "w", context)
        cls._set_create_metadata(handle)
        return handle
