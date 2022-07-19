from typing import List, Optional

# from .logging import log_io
from .tiledb_array import TileDBArray
from .tiledb_group import TileDBGroup

# from typing import Optional, Sequence, Set, Tuple

# import numpy as np
# import pandas as pd
# import tiledb

# import tiledbsc.v1.util as util


# from .types import Ids


class SOMADataFrame(TileDBArray):
    """
    Represents ``obs``, ``var``, and others.
    """

    # TODO
    # A SOMADataFrame contains a "pseudo-column" called __rowid, of type uint64 and domain
    # [0,num_rows).  The __rowid pseudo-column contains a unique value for each row in the
    # SOMADataFrame, and is intended to act as a join key for other objects, such as a SOMANdArray.

    # ----------------------------------------------------------------
    # create(string uri, Arrow.Schema schema,  user_indexed=True, string[] index_column_names) -> void
    #
    # create(string uri, Arrow.Schema schema,  user_indexed=False) -> void
    #
    # Parameters:
    #
    # * uri - location at which to create the object
    #
    # * schema -- an Arrow Schema defining the per-column schema. This schema must define all columns,
    #   including columns to be named as index columns. The column name ``__rowid`` is reserved for the
    #   pseudo-column of the same name. If the schema includes types unsupported by the SOMA
    #   implementation, an error will be raised.
    #
    # * user_indexed - boolean. If ``false``, is a ``row-indexed`` dataframe. If ``true``, is a ``user-indexed``
    #   dataframe.
    #
    # * index_column_names - a list of column names to use as user-defined index columns (e.g.,
    #   ``['cell_type', 'tissue_type']``). All named columns must exist in the schema, and at least one
    #   index column name is required. This parameter is undefined if ``user_indexed`` is False (i.e., if
    #   the dataframe is ``row-indexed``).

    _is_user_indexed: bool
    _index_column_names: List[str]

    def __init__(
        self,
        uri: str,
        name: str,
        # TODO: support row-indexed
        index_column_names: List[str],
        *,
        parent: Optional[TileDBGroup] = None,
    ):
        """
        Also see the :class:`TileDBObject` constructor.
        """
        # TODO: more options
        # assert name in ["obs", "var"]

        super().__init__(uri=uri, name=name, parent=parent)

        # TODO: add support
        self._is_user_indexed = False
        self._index_column_names = index_column_names

        # TODO DEFER self.dim_name = name + "_id"

    # TODO: static/class method?
    #    def delete(uri: str) -> None
    #        """
    #        Delete the SOMADataFrame specified with the URI.
    #        """

    # TODO: static/class method?
    #    def exists(uri: str) -> bool
    #        """
    #        Return true if object exists and is a SOMADataFrame.
    #        """

    #    def get_metadata():
    #        """
    #        Access the metadata as a mutable [`SOMAMetadataMapping`](#SOMAMetadataMapping)
    #        """

    # get_type() is inherited from TileDBObject

    #    def get_shape() -> Tuple[int]:
    #        """
    #        Return length of each dimension, always a list of length ``ndims``
    #        """

    def get_ndims(self) -> int:
        """
        Return number of index columns
        """
        return len(self._index_column_names)

    #    def get_schema(self) -> Arrow.Schema:
    #        """
    #        Return data schema, in the form of an Arrow Schema
    #        """

    def get_is_user_indexed(self) -> bool:
        """
        Return true if user-indexed, false if row-indexed.
        """
        return self._is_user_indexed

    def get_index_column_names(self) -> List[str]:
        """
        Return index (dimension) column names if user-indexed, or an empty list if row-indexed.
        """
        if self._is_user_indexed:
            return self._index_column_names
        else:
            return []


# ----------------------------------------------------------------
#    def read():
#        """
#        Read a subset of data from the SOMADataFrame
#        """

# ### Operation: read()
#
# Read a user-defined subset of data, addressed by the dataframe indexing columns, optionally
# filtered, and return results as one or more Arrow.RecordBatch.
#
# Summary:
#
# ```
# read(
#     ids=[[id,...]|all, ...],
#     column_names=[`string`, ...]|all,
#     partitions,
#     result_order,
#     value_filter
# ) -> delayed iterator over Arrow.RecordBatch
# ```
#
# Parameters:
#
# o ids - for each index dimension, which rows to read. Defaults to 'all'.
# o column_names - the named columns to read and return. Defaults to all.
# o partitions - an optional [`SOMAReadPartitions`](#SOMAReadPartitions) hint to indicate how
#   results should be organized.
# o result_order - order of read results. If dataframe is `user-indexed`, can be one of row-major,
#   col-major or unordered. If dataframe is `row-indexed`, can be one of rowid-ordered or unordered.
# o value_filter - an optional [value filter](#value-filters) to apply to the results. Defaults to
#   no filter.
#
# **Indexing**: the `ids` parameter will support per-dimension:
#
# o for `row-indexed` dataframes, a row offset (uint), a row-offset range (slice), or a list of both.
# o for `user-indexed` dataframes, a list of values of the type of the indexed column.
#
# The `read` operation will return a language-specific iterator over one or more Arrow RecordBatch
# objects, allowing the incremental processing of results larger than available memory. The actual
# iterator used is delegated to language-specific SOMA specs.

# ----------------------------------------------------------------
#    def write():
#        """
#        Write a subset of data to the SOMADataFrame
#        """

# Write an Arrow.RecordBatch to the persistent object. As duplicate index values are not allowed,
# index values already present in the object are overwritten and new index values are added.
#
# ```
# write(Arrow.RecordBatch values)
# ```
#
# Parameters:
#
# o values - an Arrow.RecordBatch containing all columns, including the index columns. The schema
#   for the values must match the schema for the SOMADataFrame.
#
# If the dataframe is `row-indexed`, the `values` Arrow RecordBatch must contain a `__rowid`
# (uint64) column, indicating which rows are being written.
