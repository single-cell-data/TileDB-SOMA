from typing import Optional, Tuple

# from .logging import log_io
from .tiledb_array import TileDBArray
from .tiledb_group import TileDBGroup

# from typing import Optional, Sequence, Set, Tuple

# import tiledb

# import tiledbsc.v1.util as util


# from .types import Ids


# TODO: rethink parenting -- add a middle layer
class SOMADenseNdArray(TileDBArray):
    """
    Represents ``X`` and others.
    """

    # ----------------------------------------------------------------
    # create(uri, ...)
    # Create a SOMADenseNdArray named with the URI.

    # ### Operation: create()
    #
    # Create a new SOMADenseNdArray with user-specified URI and schema.
    #
    # ```
    # create(string uri, type, shape) -> void
    # ```
    #
    # Parameters:
    #
    # - uri - location at which to create the object
    # - type - an Arrow type defining the type of each element in the array. If the type is unsupported, an error will be raised.
    # - shape - the length of each domain as a list, e.g., [100, 10]. All lengths must be in the uint64 range.

    _shape: Tuple[int]

    def __init__(
        self,
        uri: str,
        # TODO: incorporate
        name: str,
        # TODO: type,
        shape: Tuple[int],
        *,
        parent: Optional[TileDBGroup] = None,
    ):
        """
        Also see the :class:`TileDBObject` constructor.
        """
        # TODO: more options
        # assert name in ["obs", "var"]

        super().__init__(uri=uri, name=name, parent=parent)

        self._shape = shape

    # TODO: static/class method?
    #    def delete(uri: str) -> None
    #        """
    #        Delete the SOMADenseNdArray specified with the URI.
    #        """

    # TODO: static/class method?
    #    def exists(uri: str) -> bool
    #        """
    #        Return true if object exists and is a SOMADenseNdArray.
    #        """

    #    def get_metadata():
    #        """
    #        Access the metadata as a mutable [`SOMAMetadataMapping`](#SOMAMetadataMapping)
    #        """

    # get_type() is inherited from TileDBObject

    def get_shape(self) -> Tuple[int]:
        """
        Return length of each dimension, always a list of length ``ndims``
        """
        return self._shape

    def get_ndims(self) -> int:
        """
        Return number of index columns
        """
        return len(self._index_column_names)

    #    def get_schema(self) -> Arrow.Schema:
    #        """
    #        Return data schema, in the form of an Arrow Schema
    #        """

    def get_is_sparse(self) -> bool:
        """
        Returns ``False``.
        """
        return False


# ----------------------------------------------------------------
#    def read():
#        """
#        Read a slice of data from the SOMADenseNdArray
#        """

# ### Operation: read()
#
# Read a user-specified subset of the object, and return as one or more Arrow.Tensor.
#
# Summary:
#
# ```
# read(
#     [slice, ...],
#     partitions,
#     result_order
# ) -> delayed iterator over DenseReadResult
# ```
#
# - slice - per-dimension slice, expressed as a scalar, a range, or a list of both.
# - partitions - an optional [`SOMAReadPartitions`](#SOMAReadPartitions) hint to indicate how results should be organized.
# - result_order - order of read results. Can be one of row-major or column-major.
#
# The `read` operation will return a language-specific iterator over one or more Arrow Tensor
# objects and information describing them, allowing the incremental processing of results larger
# than available memory. The actual iterator used is delegated to language-specific SOMA specs. The
# `DenseReadResult` should include:
#
# - The coordinates of the slice (e.g., origin, shape)
# - an Arrow.Tensor with the slice values

# ----------------------------------------------------------------
#    def write():
#        """
#        Write a slice of data to the SOMADenseNdArray
#        """

# ### Operation: write()
#
# Write an Arrow.Tensor to the persistent object. As duplicate index values are not allowed, index
# values already present in the object are overwritten and new index values are added.
#
# ```
# write(coords, Arrow.Tensor values)
# ```
#
# Parameters:
#
# - coords[] - location at which to write the tensor
# - values - an Arrow.Tensor containing values to be written. The type of elements in `values` must
#   match the type of the SOMADenseNdArray.
