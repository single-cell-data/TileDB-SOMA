# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

from typing import Any

import pyarrow as pa

# This package's pybind11 code
from . import pytiledbsoma as clib
from ._managed_query import ManagedQuery
from ._soma_object import SOMAObject
from ._util import _cast_domainish, _cast_record_batch


class SOMAArray(SOMAObject):
    """Base class for all SOMAArrays: DataFrame and NDarray.

    Lifecycle:
        Maturing.
    """

    __slots__ = ()

    @property
    def schema(self) -> pa.Schema:
        """Returns data schema, in the form of an
        `Arrow Schema <https://arrow.apache.org/docs/python/generated/pyarrow.Schema.html>`_.

        Lifecycle:
            Maturing.
        """
        return self._handle.schema

    def schema_config_options(self) -> clib.PlatformSchemaConfig:
        """Returns metadata about the array schema that is not encompassed within
        the Arrow Schema, in the form of a PlatformConfig.

        Available attributes are:
            * capacity: int
            * allows_duplicates: bool
            * tile_order: str
            * cell_order: str
            * offsets_filters: str
                * name (of filter): str
                * compression_level: str
            * validity_filters: str
            * attrs: str
                * name (of attribute): str
                    * filters: str
                        * name (of filter): str
                        * compression_level: str
            * dims: str
                * name (of dimension): str
                    * filters: str
                        * name (of filter): str
                        * compression_level: str
                    * tile: int

        Lifecycle:
            Experimental.
        """
        return self._handle.schema_config_options()

    def non_empty_domain(self) -> tuple[tuple[Any, Any], ...]:
        """Retrieves the non-empty domain for each dimension, namely the smallest
        and largest indices in each dimension for which the array/dataframe has
        data occupied.  This is nominally the same as the domain used at
        creation time, but if for example only a portion of the available domain
        has actually had data written, this function will return a tighter
        range.
        """
        return _cast_domainish(self._handle.non_empty_domain())

    def _tiledb_array_keys(self) -> tuple[str, ...]:
        """Return all dim and attr names."""
        return self._tiledb_dim_names() + self._tiledb_attr_names()

    def _tiledb_dim_names(self) -> tuple[str, ...]:
        """Reads the dimension names from the schema: for example, ['obs_id', 'var_id']."""
        return tuple(self._handle.dimension_names)

    def _tiledb_attr_names(self) -> tuple[str, ...]:
        """Reads the attribute names from the schema:
        for example, the list of column names in a dataframe.
        """
        return tuple(f.name for f in self.schema if f.name not in self._handle.dimension_names)

    def _domain(self) -> tuple[tuple[Any, Any], ...]:
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
        return _cast_domainish(self._handle.domain())

    def _maxdomain(self) -> tuple[tuple[Any, Any], ...]:
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
        return _cast_domainish(self._handle.maxdomain())

    def _write_table(self, values: pa.Table, sort_coords: bool) -> None:
        """Helper function that sets the correct result order for the layout
        and allows for multiple submissions before calling `submit` and `finalize`
        for unordered write or `submit_and_finalize` for global order writes.

        Args:
            values:
                An `Arrow table <https://arrow.apache.org/docs/python/generated/pyarrow.Table.html>`_
                containing all columns, including the index columns. The schema
                for the values must match the schema for the :class:`DataFrame`.

                If a column is of categorical type in the schema and a
                flattened/non-categorical column is presented for data on write,
                a ``ValueError`` is raised.  If a column is of non-categorical
                type in the schema and a categorical column is presented for data
                on write, the data are written as an array of category values,
                and the category-type information is not saved.
            sort_coords:
                Whether the coordinates need to be sorted (True) or are already
                sorted in global order (False). In the PlatformConfig, this is
                is to True by default.
        """
        batches = values.to_batches()
        if not batches:
            return

        array_schema = self.schema
        for name in values.schema.names:
            if name not in array_schema.names:
                raise ValueError(
                    f"Cannot write data. Field '{name}' in the input data is not a column in this {self._handle_type.__name__}."
                )
        batch_schema = pa.schema([array_schema.field(name) for name in values.schema.names])

        if sort_coords:
            # Finalize each batch as it is written.
            for batch in batches:
                mq = ManagedQuery(self)
                mq._handle.set_layout(clib.ResultOrder.unordered)
                mq.submit_batch(_cast_record_batch(batch, batch_schema, safe=True))
                mq.finalize()
        else:
            # Single global order query - only finalize at the end.
            mq = ManagedQuery(self)
            mq._handle.set_layout(clib.ResultOrder.globalorder)
            for batch in batches[:-1]:
                mq.submit_batch(_cast_record_batch(batch, batch_schema, safe=True))
            mq._handle.submit_and_finalize_batch(_cast_record_batch(batches[-1], batch_schema, safe=True))
