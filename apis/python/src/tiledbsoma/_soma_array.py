# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

import warnings
from typing import Any, Tuple

import pyarrow as pa

from . import _tdb_handles

# This package's pybind11 code
from . import pytiledbsoma as clib  # noqa: E402
from ._read_iters import ManagedQuery
from ._soma_object import SOMAObject


class SOMAArray(SOMAObject[_tdb_handles.SOMAArrayWrapper[Any]]):
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

    def config_options_from_schema(self) -> clib.PlatformConfig:
        """Returns metadata about the array that is not encompassed within the
        Arrow Schema, in the form of a PlatformConfig (deprecated).

        Use ``schema_config_options`` instead.

        Available attributes are:
            * dataframe_dim_zstd_level: int
            * sparse_nd_array_dim_zstd_level: int
            * sparse_nd_array_dim_zstd_level: int
            * write_X_chunked: bool
            * goal_chunk_nnz: int
            * remote_cap_nbytes: int
            * capacity: int
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
            * allows_duplicates: bool
            * tile_order: str
            * cell_order: str
            * consolidate_and_vacuum: bool

        Lifecycle:
            Deprecated.
        """
        warnings.warn(
            "Deprecated. Use schema_config_options instead.", DeprecationWarning
        )
        return self._handle.config_options_from_schema()

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

        layout = (
            clib.ResultOrder.unordered if sort_coords else clib.ResultOrder.globalorder
        )

        if layout == clib.ResultOrder.unordered:
            for batch in batches:
                # Create new ManagedQuery per each batch
                mq = ManagedQuery(self)._handle
                mq.set_layout(layout)
                mq.set_array_data(batch)

                # Submit and flush every batch
                mq.submit_write()
                mq.finalize()

        else:  # global order
            # Create a single ManagedQuery at the beginning
            mq = ManagedQuery(self)._handle
            mq.set_layout(layout)

            # Submit for each batch but don't finalize
            for batch in batches[:-1]:
                mq.set_array_data(batch)
                mq.submit_write()

            # Only finalize at the last batch
            mq.set_array_data(batches[-1])
            mq.submit_and_finalize()
