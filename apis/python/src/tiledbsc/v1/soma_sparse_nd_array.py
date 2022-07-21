from typing import Optional, Tuple

# from .logging import log_io
from .tiledb_array import TileDBArray
from .tiledb_group import TileDBGroup

# from typing import Optional, Sequence, Set, Tuple

# import tiledb

# import tiledbsc.v1.util as util

# TODO: deconflate the CREATE method and the CONSTRUCTOR.
# These are not the same.


class SOMASparseNdArray(TileDBArray):
    """
    Represents ``X`` and others.
    """

    def __init__(
        self,
        uri: str,
        *,
        name: Optional[str] = None,
        parent: Optional[TileDBGroup] = None,
    ):
        """
        Also see the :class:`TileDBObject` constructor.
        """

        super().__init__(uri=uri, name=name, parent=parent)

    # ```
    # create(string uri, type, shape) -> void
    # ```
    #
    # Parameters:
    #
    # - uri - location at which to create the object
    #
    # - type - an Arrow type defining the type of each element in the array. If the type is
    #   unsupported, an error will be raised.
    #
    # - shape - the length of each domain as a list, e.g., [100, 10]. All lengths must be in the
    #   uint64 range.
    #        Create a new SOMASparseNdArray with user-specified URI and schema.

    #        # ----------------------------------------------------------------
    #        # TODO: type,
    #        shape: Tuple,
    #
    #        # Check that ndims, and each dimension, are positive
    #        assert len(shape) > 0
    #        for e in shape:
    #            assert e > 0

    # TODO: static/class method?
    #    def delete(uri: str) -> None
    #        """
    #        Delete the SOMASparseNdArray specified with the URI.
    #        """

    # TODO: static/class method?
    #    def exists(uri: str) -> bool
    #        """
    #        Return true if object exists and is a SOMASparseNdArray.
    #        """

    #    def get_metadata():
    #        """
    #        Access the metadata as a mutable [`SOMAMetadataMapping`](#SOMAMetadataMapping)
    #        """

    # get_type() is inherited from TileDBObject

    def get_shape(self) -> Tuple:
        """
        Return length of each dimension, always a list of length ``ndims``
        """
        return self._shape

    def get_ndims(self) -> int:
        """
        Return number of index columns
        """
        return len(self._shape)

    #    def get_schema(self) -> Arrow.Schema:
    #        """
    #        Return data schema, in the form of an Arrow Schema
    #        """

    def get_is_sparse(self) -> bool:
        """
        Returns ``True``.
        """
        return True


#    def get_nnz(self) -> wint:
#        """
#        Return the number of non-zero values in the array
#        """
#        return 999


# ----------------------------------------------------------------
#    def read():
#        """
#        Read a slice of data from the SOMASparseNdArray
#        """

# ### Operation: read()
#
# Read a user-specified subset of the object, and return as one or more Arrow.SparseTensor.
#
# Summary:
#
# ```
# read(
#     [slice, ...],
#     partitions,
#     result_order
# ) -> delayed iterator over Arrow.SparseTensor
# ```
#
# - slice - per-dimension slice, expressed as a scalar, a range, or a list of both.
# - partitions - an optional [`SOMAReadPartitions`](#SOMAReadPartitions) hint to indicate how
#   results should be organized.
# - result_order - order of read results. Can be one of row-major, column-major and unordered.
#
# The `read` operation will return a language-specific iterator over one or more Arrow SparseTensor
# objects, allowing the incremental processing of results larger than available memory. The actual
# iterator used is delegated to language-specific SOMA specs.

# ----------------------------------------------------------------
#    def write():
#        """
#        Write a slice of data to the SOMASparseNdArray
#        """

# ### Operation: write()
#
# ```
# write(Arrow.SparseTensor values)
# ```
#
# Parameters:
#
# - values - an Arrow.SparseTensor containing values to be written. The type of elements in `values`
#   must match the type of the SOMASparseNdArray.
