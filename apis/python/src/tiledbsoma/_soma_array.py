# Copyright (c) 2021-2023 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2023 TileDB, Inc.
#
# Licensed under the MIT License.

from typing import Any, Optional, Sequence, Tuple

import pyarrow as pa
from somacore import options
from typing_extensions import Self

from . import _tdb_handles, _util

# This package's pybind11 code
from . import pytiledbsoma as clib  # noqa: E402
from ._soma_object import SOMAObject
from ._types import OpenTimestamp, is_nonstringy_sequence
from .options._soma_tiledb_context import SOMATileDBContext


class SOMAArray(SOMAObject[_tdb_handles.SOMAArrayWrapper[Any]]):
    """Base class for all SOMAArrays: DataFrame and NDarray.

    Lifecycle:
        Maturing.
    """

    __slots__ = ()

    @classmethod
    def open(
        cls,
        uri: str,
        mode: options.OpenMode = "r",
        *,
        tiledb_timestamp: Optional[OpenTimestamp] = None,
        context: Optional[SOMATileDBContext] = None,
        platform_config: Optional[options.PlatformConfig] = None,
        clib_type: Optional[str] = None,
    ) -> Self:
        """Opens this specific type of SOMA object."""
        return super().open(
            uri,
            mode,
            tiledb_timestamp=tiledb_timestamp,
            context=context,
            platform_config=platform_config,
            clib_type="SOMAArray",
        )

    @property
    def schema(self) -> pa.Schema:
        """Returns data schema, in the form of an
        `Arrow Schema <https://arrow.apache.org/docs/python/generated/pyarrow.Schema.html>`_.

        Lifecycle:
            Maturing.
        """
        return self._handle.schema

    def non_empty_domain(self) -> Tuple[Tuple[Any, Any], ...]:
        """
        Retrieves the non-empty domain for each dimension, namely the smallest
        and largest indices in each dimension for which the array/dataframe has
        data occupied.  This is nominally the same as the domain used at
        creation time, but if for example only a portion of the available domain
        has actually had data written, this function will return a tighter
        range.
        """
        return self._handle.non_empty_domain()

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

    def _domain(self) -> Tuple[Tuple[Any, Any], ...]:
        """This is the SOMA domain, not the core domain.
        * For arrays with core current-domain support:
          o soma domain is core current domain
          o soma maxdomain is core domain
        * For arrays without core current-domain support:
          o soma domain is core domain
          o soma maxdomain is core domain
          o core current domain is not accessed at the soma level
        * Core domain has been around forever and is immutable.
        * Core current domain is new as of core 2.25 and can be
          resized up to core (max) domain.
        """
        return self._handle.domain

    def _maxdomain(self) -> Tuple[Tuple[Any, Any], ...]:
        """This is the SOMA maxdomain, not the core domain.
        * For arrays with core current-domain support:
          o soma domain is core current domain
          o soma maxdomain is core domain
        * For arrays without core current-domain support:
          o soma domain is core domain
          o soma maxdomain is core domain
          o core current domain is not accessed at the soma level
        * Core domain has been around forever and is immutable.
        * Core current domain is new as of core 2.25 and can be
          resized up to core (max) domain.
        """
        return self._handle.maxdomain

    def _set_reader_coords(self, sr: clib.SOMAArray, coords: Sequence[object]) -> None:
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
        self, sr: clib.SOMAArray, dim_idx: int, dim: pa.Field, coord: object
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
