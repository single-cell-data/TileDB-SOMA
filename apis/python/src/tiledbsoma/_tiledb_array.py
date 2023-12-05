# Copyright (c) 2021-2023 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2023 TileDB, Inc.
#
# Licensed under the MIT License.

import ctypes
import os
import sys
from typing import Any, List, Optional, Sequence, Tuple

import pyarrow as pa
import tiledb

from . import _tdb_handles, _util
from ._arrow_types import tiledb_schema_to_arrow
from ._tdb_handles import SOMAArray
from ._tiledb_object import TileDBObject
from ._types import OpenTimestamp, is_nonstringy_sequence
from .options._soma_tiledb_context import SOMATileDBContext


def _load_libs() -> None:
    """Loads the required TileDB-SOMA native library."""
    if sys.platform == "darwin":
        lib_name = "libtiledbsoma.dylib"
    else:
        lib_name = "libtiledbsoma.so"

    try:
        # Try loading the bundled native library.
        lib_dir = os.path.dirname(os.path.abspath(__file__))
        ctypes.CDLL(os.path.join(lib_dir, lib_name))
    except OSError:
        # Otherwise try loading by name only.
        ctypes.CDLL(lib_name)


# Load native libraries
_load_libs()

# This package's pybind11 code


class TileDBArray(TileDBObject[_tdb_handles.ArrayWrapper]):
    """Wraps arrays from TileDB-Py by retaining a URI, options, etc.
    Also serves as an abstraction layer to hide TileDB-specific details
    from the API, unless requested.

    Lifecycle:
        Experimental.
    """

    __slots__ = ()

    _wrapper_type = _tdb_handles.ArrayWrapper

    @property
    def schema(self) -> pa.Schema:
        """Returns data schema, in the form of an
        `Arrow Schema <https://arrow.apache.org/docs/python/generated/pyarrow.Schema.html>`_.

        Lifecycle:
            Experimental.
        """
        if isinstance(self._tiledb_array_schema(), tiledb.ArraySchema):
            return tiledb_schema_to_arrow(
                self._tiledb_array_schema(), self.uri, self._ctx
            )
        else:
            return self._tiledb_array_schema()

    def non_empty_domain(self) -> Optional[Tuple[Tuple[Any, Any], ...]]:
        """
        Retrieves the non-empty domain for each dimension, namely the smallest
        and largest indices in each dimension for which the array/dataframe has
        data occupied.  This is nominally the same as the domain used at
        creation time, but if for example only a portion of the available domain
        has actually had data written, this function will return a tighter
        range.
        """
        return self._handle.non_empty_domain()

    def _tiledb_array_schema(self) -> tiledb.ArraySchema:
        """Returns the TileDB array schema, for internal use."""
        return self._handle.schema

    def _tiledb_array_keys(self) -> Tuple[str, ...]:
        """Return all dim and attr names."""
        return self._tiledb_dim_names() + self._tiledb_attr_names()

    def _tiledb_dim_names(self) -> Tuple[str, ...]:
        """Reads the dimension names from the schema: for example, ['obs_id', 'var_id']."""
        return self._handle.dim_names

    def _tiledb_attr_names(self) -> Tuple[str, ...]:
        """Reads the attribute names from the schema:
        for example, the list of column names in a dataframe.
        """
        return self._handle.attr_names

    def _tiledb_domain(self) -> Tuple[Tuple[Any, Any], ...]:
        return self._handle.domain

    def _set_reader_coords(self, sr: SOMAArray, coords: Sequence[object]) -> None:
        """Parses the given coords and sets them on the SOMA Reader."""
        if not is_nonstringy_sequence(coords):
            raise TypeError(
                f"coords type {type(coords)} must be a regular sequence,"
                " not str or bytes"
            )

        if len(coords) > self._handle.ndim:
            raise ValueError(
                f"coords ({len(coords)} elements) must be shorter than ndim"
                f" ({self._handle.ndim})"
            )
        for i, coord in enumerate(coords):
            dim = self.schema.field(i)
            if not self._set_reader_coord(sr, i, dim, coord):
                raise TypeError(
                    f"coord type {type(coord)} for dimension {dim.name}"
                    f" (slot {i}) unsupported"
                )

    def _set_reader_coord(
        self, sr: SOMAArray, dim_idx: int, dim: pa.Field, coord: object
    ) -> bool:
        """Parses a single coordinate entry.

        The base implementation parses the most fundamental types shared by all
        TileDB Array types; subclasses can implement their own readers that
        handle types not recognized here.

        Returns:
            True if successful, False if unrecognized.
        """
        if coord is None:
            return True  # No constraint; select all in this dimension

        if isinstance(coord, int):
            sr.set_dim_points_int64(dim.name, [coord])
            return True
        if isinstance(coord, slice):
            _util.validate_slice(coord)
            try:
                dom = self._handle.domain[dim_idx]
                lo_hi = _util.slice_to_numeric_range(coord, dom)
            except _util.NonNumericDimensionError:
                return False  # We only handle numeric dimensions here.
            if lo_hi:
                sr.set_dim_ranges_int64(dim.name, [lo_hi])
            # If `None`, coord was `slice(None)` and there is no constraint.
            return True
        return False

    @classmethod
    def _create_internal(
        cls,
        uri: str,
        schema: tiledb.ArraySchema,
        context: SOMATileDBContext,
        tiledb_timestamp: Optional[OpenTimestamp],
    ) -> _tdb_handles.ArrayWrapper:
        """Creates the TileDB Array for this type and returns an opened handle.

        This does the work of creating a TileDB Array with the provided schema
        at the given URI, sets the necessary metadata, and returns a handle to
        the newly-created array, open for writing.
        """
        tiledb.Array.create(uri, schema, ctx=context.tiledb_ctx)
        handle = cls._wrapper_type.open(uri, "w", context, tiledb_timestamp)
        cls._set_create_metadata(handle)
        return handle

    def _consolidate_and_vacuum(
        self, modes: List[str] = ["fragment_meta", "commits"]
    ) -> None:
        """
        This post-ingestion helper consolidates and vacuums fragment metadata and commit files --
        this is quick to do, and positively impacts query performance.  It does _not_ consolidate
        bulk array data, which is more time-consuming and should be done at the user's opt-in
        discretion.
        """

        for mode in modes:
            self._consolidate(modes=[mode])
            self._vacuum(modes=[mode])

    def _consolidate(self, modes: List[str] = ["fragment_meta", "commits"]) -> None:
        """
        This post-ingestion helper consolidates by default fragment metadata and commit files --
        this is quick to do, and positively impacts query performance.
        """

        for mode in modes:
            cfg = self._ctx.config()
            cfg["sm.consolidation.mode"] = mode
            ctx = tiledb.Ctx(cfg)

            tiledb.consolidate(self.uri, ctx=ctx)

    def _vacuum(self, modes: List[str] = ["fragment_meta", "commits"]) -> None:
        """
        This post-ingestion helper vacuums by default fragment metadata and commit files. Vacuuming is not multi-process safe and requires coordination that nothing is currently reading the files that will be vacuumed.
        """

        for mode in modes:
            cfg = self._ctx.config()
            cfg["sm.vacuum.mode"] = mode
            ctx = tiledb.Ctx(cfg)

            tiledb.vacuum(self.uri, ctx=ctx)
