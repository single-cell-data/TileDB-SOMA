# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.
"""Definitions of data storage interfaces for SOMA implementations.

SOMA users should ordinarily not need to import this module directly; relevant
members will be exported to the ``somacore`` namespace.

Default values are provided here as a reference for implementors.
"""

from __future__ import annotations

import abc
from collections.abc import Iterator, Sequence
from typing import (
    Any,
    ClassVar,
    Final,
    TypeVar,
    Union,
)

import pyarrow as pa
from typing_extensions import Self

from . import base, options
from .types import StatusAndReason

_RO_AUTO = options.ResultOrder.AUTO

AxisDomain = Union[tuple[Any, Any], list[Any], None]
Domain = Sequence[AxisDomain]

_UNBATCHED = options.BatchSize()


class DataFrame(base.SOMAObject, metaclass=abc.ABCMeta):
    """A multi-column table with a user-defined schema.

    Lifecycle: maturing
    """

    __slots__ = ()
    soma_type: Final = "SOMADataFrame"  # type: ignore[misc]

    # Lifecycle

    @classmethod
    @abc.abstractmethod
    def create(
        cls,
        uri: str,
        *,
        schema: pa.Schema,
        index_column_names: Sequence[str] = (options.SOMA_JOINID,),
        domain: Sequence[tuple[Any, Any] | None] | None = None,
        platform_config: options.PlatformConfig | None = None,
        context: Any | None = None,  # noqa: ANN401
    ) -> Self:
        """Creates a new ``DataFrame`` at the given URI.

        The schema of the created dataframe will include a column named
        ``soma_joinid`` of type ``pyarrow.int64``, with negative values
        disallowed.  If a ``soma_joinid`` column is present in the provided
        schema, it must be of the correct type.  If no ``soma_joinid`` column
        is provided, one will be added.  It may be used as an indexed column.

        Args:
            uri: The URI where the dataframe will be created.

            schema: Arrow schema defining the per-column schema. This schema
                must define all columns, including columns to be named as index
                columns.  If the schema includes types unsupported by the SOMA
                implementation, an error will be raised.

            index_column_names: A list of column names to use as user-defined
                index columns (e.g., ``['cell_type', 'tissue_type']``).
                All named columns must exist in the schema, and at least one
                index column name is required.

            domain:
                An optional sequence of tuples specifying the domain of each
                index column. Each tuple must be a pair consisting of the
                minimum and maximum values storable in the index column. For
                example, if there is a single int64-valued index column, then
                ``domain`` might be ``[(100, 200)]`` to indicate that values
                between 100 and 200, inclusive, can be stored in that column.
                If provided, this sequence must have the same length as
                ``index_column_names``, and the index-column domain will be as
                specified.  If omitted entirely, or if ``None`` in a given
                dimension, the corresponding index-column domain will use an
                empty range, and data writes after that will fail with an
                exception.  Unless you have a particular reason not to, you
                should always provide the desired `domain` at create time: this
                is an optional but strongly recommended parameter. See also
                ``change_domain`` which allows you to expand the domain after
                create.

            platform_config: platform-specific configuration; keys are SOMA
                implementation names.

            context: Other implementation-specific configuration.

        Returns:
            The newly created dataframe, opened for writing.

        Lifecycle: maturing
        """
        raise NotImplementedError

    # Data operations

    @abc.abstractmethod
    def read(
        self,
        coords: options.SparseDFCoords = (),
        column_names: Sequence[str] | None = None,
        *,
        batch_size: options.BatchSize = _UNBATCHED,
        partitions: options.ReadPartitions | None = None,
        result_order: options.ResultOrderStr = _RO_AUTO,
        value_filter: str | None = None,
        platform_config: options.PlatformConfig | None = None,
    ) -> ReadIter[pa.Table]:
        """Reads a user-defined slice of data into Arrow tables.

        Args:
            coords: for each index dimension, which rows to read.
                Defaults to ``()``, meaning no constraint -- all IDs.
            column_names: the named columns to read and return.
                Defaults to ``None``, meaning no constraint -- all column names.
            partitions: If present, specifies that this is part of
                a partitioned read, and which part of the data to include.
            result_order: the order to return results, specified as a
                :class:`~options.ResultOrder` or its string value.
            value_filter: an optional value filter to apply to the results.
                The default of ``None`` represents no filter. Value filter
                syntax is implementation-defined; see the documentation
                for the particular SOMA implementation for details.

        Returns:
            A :class:`ReadIter` of :class:`pa.Table`s.


        **Indexing:**

        Indexing is performed on a per-column basis for each indexed column.
        To specify dimensions:

        - A sequence of coordinates is accepted, one per indexed dimension.
        - The sequence length must be less than or equal to the number of
          indexed dimensions.
        - If the sequence is shorter than the number of indexed coordinates,
          then no constraint (i.e. ``None``) is used for the remaining
          indexed dimensions.
        - Specifying an empty sequence (e.g. ``()``, the default) represents
          no constraints over any dimension, returning the entire dataset.

        Each dimension may be indexed as follows:

        - ``None`` or ``slice(None)`` places no constraint on the dimension.
        - Coordinates can be specified as a scalar value, a Python sequence
          (``list``, ``tuple``, etc.), a NumPy ndarray, an Arrow array, or
          similar objects (as defined by ``SparseDFCoords``).
        - Slices specify a closed range: ``slice(2, 4)`` includes both 2 and 4.
          Slice *steps* may not be used: ``slice(10, 20, 2)`` is invalid.
          ``slice(None)`` places no constraint on the dimension. Half-specified
          slices like ``slice(None, 99)`` and ``slice(5, None)`` specify
          all indices up to and including the value, and all indices
          starting from and including the value.
        - Negative values in indices and slices are treated as raw domain values
          and not as indices relative to the end, unlike traditional Python
          sequence indexing.

        Lifecycle: maturing
        """
        raise NotImplementedError

    @abc.abstractmethod
    def change_domain(
        self,
        newdomain: Domain,
        check_only: bool = False,
    ) -> StatusAndReason:
        """Allows you to enlarge the domain of a SOMA :class:`DataFrame`, when
        the ``DataFrame`` already has a domain.

        The argument must be a tuple of pairs of low/high values for the desired
        domain, one pair per index column. For string index columns, you must
        offer the low/high pair as `("", "")`, or as ``None``.  If ``check_only``
        is ``True``, returns whether the operation would succeed if attempted,
        and a reason why it would not.

        For example, suppose the dataframe's sole index-column name is
        ``"soma_joinid"`` (which is the default at create).  If the dataframe's
        ``.maxdomain`` is ``((0, 999999),)`` and its ``.domain`` is ``((0,
        2899),)``, this means that ``soma_joinid`` values between 0 and 2899 can
        be read or written; any attempt to read or write ``soma_joinid`` values
        outside this range will result in an error. If you then apply
        ``.change_domain([(0, 5700)])``, then ``.domain`` will
        report ``((0, 5699),)``, and now ``soma_joinid`` values in the range 0
        to 5699 can now be written to the dataframe.

        If you use non-default ``index_column_names`` in the dataframe's
        ``create`` then you need to specify the (low, high) pairs for each
        index column. For example, if the dataframe's ``index_column_names``
        is ``["soma_joinid", "cell_type"]``, then you can upgrade domain using
        ``[(0, 5699), ("", "")]``.

        Lastly, it is an error to try to set the ``domain`` to be smaller than
        ``maxdomain`` along any index column.  The ``maxdomain`` of a dataframe is
        set at creation time, and cannot be extended afterward.

        Lirecycle: maturing
        """
        raise NotImplementedError

    @abc.abstractmethod
    def write(
        self,
        values: pa.RecordBatch | pa.Table,
        *,
        platform_config: options.PlatformConfig | None = None,
    ) -> Self:
        """Writes the data from an Arrow table to the persistent object.

        As duplicate index values are not allowed, index values already present
        in the object are overwritten and new index values are added.

        Args:
            values: An Arrow table containing all columns, including
                the index columns. The schema for the values must match
                the schema for the ``DataFrame``.

        Returns: ``self``, to enable method chaining.

        Lifecycle: maturing
        """
        raise NotImplementedError

    # Metadata operations

    @property
    @abc.abstractmethod
    def schema(self) -> pa.Schema:
        """The schema of the data in this dataframe.

        Lifecycle: maturing
        """
        raise NotImplementedError

    @property
    @abc.abstractmethod
    def index_column_names(self) -> tuple[str, ...]:
        """The names of the index (dimension) columns.

        Lifecycle: maturing
        """
        raise NotImplementedError

    @property
    @abc.abstractmethod
    def domain(self) -> tuple[tuple[Any, Any], ...]:
        """The allowable range of values in each index column.

        Returns: a tuple of minimum and maximum values, inclusive,
            storable on each index column of the dataframe.

        Lifecycle: maturing
        """
        raise NotImplementedError


class NDArray(base.SOMAObject, metaclass=abc.ABCMeta):
    """Common behaviors of N-dimensional arrays of a single primitive type."""

    __slots__ = ()

    # Lifecycle

    @classmethod
    @abc.abstractmethod
    def create(
        cls,
        uri: str,
        *,
        type: pa.DataType,
        shape: Sequence[int | None],
        platform_config: options.PlatformConfig | None = None,
        context: Any | None = None,  # noqa: ANN401
    ) -> Self:
        """Creates a new ND array of the current type at the given URI.

        Args:
            uri: The URI where the array will be created.
            type: The Arrow type to store in the array.
                If the type is unsupported, an error will be raised.
            shape: The maximum capacity of each dimension, including room
                for any intended future appends, specified as one element
                per dimension, e.g. ``(100, 10)``.  All lengths must be in
                the positive int64 range, or ``None``.  It's necessary to say
                ``shape=(None, None)`` or ``shape=(None, None, None)``,
                as the sequence length determines the number of dimensions
                (N) in the N-dimensional array.

                For sparse arrays only, if a slot is None, then the maximum
                possible int64 will be used, making a sparse array growable.

        Returns: The newly created array, opened for writing.

        Lifecycle: maturing
        """
        raise NotImplementedError

    def resize(self, newshape: Sequence[int | None], check_only: bool = False) -> StatusAndReason:
        """Increases the shape of the array as specfied. Raises an error if the new
        shape is less than the current shape in any dimension. Raises an error if
        the new shape exceeds maxshape in any dimension. Raises an error if the
        array doesn't already have a shape: in that case please call
        tiledbsoma_upgrade_shape. If ``check_only`` is ``True``, returns
        whether the operation would succeed if attempted, and a reason why it
        would not.

        Lifecycle: maturing
        """
        raise NotImplementedError

    # Metadata operations

    @property
    @abc.abstractmethod
    def shape(self) -> tuple[int, ...]:
        """The maximum capacity (domain) of each dimension of this array.

        Lifecycle: maturing
        """
        raise NotImplementedError

    @property
    def ndim(self) -> int:
        """The number of dimensions in this array.

        Lifecycle: maturing
        """
        return len(self.shape)

    @property
    @abc.abstractmethod
    def schema(self) -> pa.Schema:
        """The schema of the data in this array.

        Lifecycle: maturing
        """
        raise NotImplementedError

    is_sparse: ClassVar[bool]
    """True if the array is sparse, False if it is dense.

    Lifecycle: maturing
    """


class DenseNDArray(NDArray, metaclass=abc.ABCMeta):
    """An N-dimensional array stored densely.

    Lifecycle: maturing
    """

    __slots__ = ()
    soma_type: Final = "SOMADenseNDArray"  # type: ignore[misc]
    is_sparse: Final = False  # type: ignore[misc]

    @abc.abstractmethod
    def read(
        self,
        coords: options.DenseNDCoords = (),
        *,
        partitions: options.ReadPartitions | None = None,
        result_order: options.ResultOrderStr = _RO_AUTO,
        platform_config: options.PlatformConfig | None = None,
    ) -> pa.Tensor:
        """Reads the specified subarray as a Tensor.

        Coordinates must specify a contiguous subarray, and the number of
        coordinates must be less than or equal to the number of dimensions.

        Args:
            coords: A per-dimension sequence of coordinates defining
                the range to read.
            partitions: If present, specifies that this is part of
                a partitioned read, and which part of the data to include.
            result_order: the order to return results, specified as a
                :class:`~options.ResultOrder` or its string value.
            platform_config: platform-specific configuration; keys are SOMA
                implementation names.

        Returns: The data over the requested range as a tensor.

        **Indexing:**

        Indexing is performed on a per-dimension basis.

        - A sequence of coordinates is accepted, one per dimension.
        - The sequence length must be less than or equal to
          the number of dimensions.
        - If the sequence is shorter than the number of dimensions, the
          remaining dimensions are unconstrained.
        - Specifying an empty sequence (e.g. ``()``, the default) represents
          no constraints over any dimension, returning the entire dataset.

        Each dimension may be indexed by value or slice:

        - Slices specify a closed range: ``slice(2, 4)`` includes 2, 3, and 4.
          Slice *steps* may not be used: ``slice(10, 20, 2)`` is invalid.
          ``slice(None)`` places no constraint on the dimension. Half-specified
          slices like ``slice(None, 99)`` and ``slice(5, None)`` specify
          all indices up to and including the value, and all indices
          starting from and including the value.
        - Negative indexing is not supported.

        Lifecycle: maturing
        """
        raise NotImplementedError

    @abc.abstractmethod
    def write(
        self,
        coords: options.DenseNDCoords,
        values: pa.Tensor,
        *,
        platform_config: options.PlatformConfig | None = None,
    ) -> Self:
        """Writes an Arrow tensor to a subarray of the persistent object.

        The subarray written is defined by ``coords`` and ``values``. This will
        overwrite existing values in the array.

        Args:
            coords: A per-dimension tuple of scalars or slices
                defining the bounds of the subarray to be written.
                See :meth:`read` for details about indexing.
            values: The values to be written to the subarray.  Must have
                the same shape as ``coords``, and matching type to the array.
        Returns: ``self``, to enable method chaining.
        Lifecycle: maturing
        """
        raise NotImplementedError


SparseArrowData = Union[
    pa.SparseCSCMatrix,
    pa.SparseCSRMatrix,
    pa.SparseCOOTensor,
    pa.Table,
]
"""Any of the sparse data storages provided by Arrow."""


class SparseNDArray(NDArray, metaclass=abc.ABCMeta):
    """A N-dimensional array stored sparsely.

    Lifecycle: maturing
    """

    __slots__ = ()
    soma_type: Final = "SOMASparseNDArray"  # type: ignore[misc]
    is_sparse: Final = True  # type: ignore[misc]

    @abc.abstractmethod
    def read(
        self,
        coords: options.SparseNDCoords = (),
        *,
        batch_size: options.BatchSize = _UNBATCHED,
        partitions: options.ReadPartitions | None = None,
        result_order: options.ResultOrderStr = _RO_AUTO,
        platform_config: options.PlatformConfig | None = None,
    ) -> SparseRead:
        """Reads the specified subarray in batches.

        Values returned are a :class:`SparseRead` object which can be converted
        to any number of formats::

            some_dense_array.read(...).tables()
            # -> an iterator of Arrow Tables

        Args:
            coords: A per-dimension sequence of coordinates defining
                the range to be read.
            batch_size: The size of batches that should be returned from a read.
            See :class:`options.BatchSize` for details.
            partitions: Specifies that this is part of a partitioned read,
                and which partition to include, if present.
            result_order: the order to return results, specified as a
                :class:`~options.ResultOrder` or its string value.
            platform_config: platform-specific configuration; keys are SOMA
                implementation names.

        Returns: The data that was requested in a :class:`SparseRead`,
            allowing access in any supported format.

        **Indexing:**

        Indexing is performed on a per-dimension basis.

        - A sequence of coordinates is accepted, one per dimension.
        - The sequence length must be less than or equal to
          the number of dimensions.
        - If the sequence is shorter than the number of dimensions, the
          remaining dimensions are unconstrained.
        - Specifying an empty sequence (e.g. ``()``, the default) represents
          no constraints over any dimension, returning the entire dataset.

        Each dimension may be indexed as follows:

        - ``None`` or ``slice(None)`` places no constraint on the dimension.
        - Coordinates can be specified as a scalar value, a Python sequence
          (``list``, ``tuple``, etc.), a ``ndarray``, an Arrow array, and
          similar objects (as defined by ``SparseNDCoords``).
        - Slices specify a closed range: ``slice(2, 4)`` includes 2, 3, and 4.
          Slice *steps* may not be used: ``slice(10, 20, 2)`` is invalid.
          ``slice(None)`` places no constraint on the dimension. Half-specified
          slices like ``slice(None, 99)`` and ``slice(5, None)`` specify
          all indices up to and including the value, and all indices
          starting from and including the value.
        - Negative indexing is not supported.

        Lifecycle: maturing
        """

    @abc.abstractmethod
    def write(
        self,
        values: SparseArrowData,
        *,
        platform_config: options.PlatformConfig | None = None,
    ) -> Self:
        """Writes a Tensor to a subarray of the persistent object.

        Args:
            values: The values to write to the array. Supported types are:

                Arrow sparse tensor: the coordinates in the tensor are
                interpreted as the coordinates to write to.  Supports the
                *experimental* types SparseCOOTensor, SparseCSRMatrix, and
                SparseCSCMatrix. There is currently no support for
                SparseCSFTensor or dense Tensor.

                Arrow table: a COO table, with columns named ``soma_dim_0``,
                    ..., ``soma_dim_N`` and ``soma_data``.
            platform_config: platform-specific configuration; keys are SOMA
                implementation names.

        Returns: ``self``, to enable method chaining.

        Lifecycle: maturing
        """
        raise NotImplementedError

    @property
    def nnz(self) -> int:
        """The number of values stored in the array, including explicit zeros.

        Lifecycle: maturing
        """
        raise NotImplementedError


#
# Read types
#

_T = TypeVar("_T")


# Sparse reads are returned as an iterable structure:


class ReadIter(Iterator[_T], metaclass=abc.ABCMeta):
    """SparseRead result iterator allowing users to flatten the iteration.

    Lifecycle: maturing
    """

    __slots__ = ()

    # __iter__ is already implemented as `return self` in Iterator.
    # SOMA implementations must implement __next__.

    @abc.abstractmethod
    def concat(self) -> _T:
        """Returns all the requested data in a single operation.

        If some data has already been retrieved using ``next``, this will return
        the remaining data, excluding that which as already been returned.

        Lifecycle: maturing
        """
        raise NotImplementedError


class SparseRead:
    """Intermediate type to choose result format when reading a sparse array.

    A query may not be able to return all of these formats. The concrete result
    may raise a ``NotImplementedError`` or may choose to raise a different
    exception (likely a ``TypeError``) containing more specific information
    about why the given format is not supported.

    Lifecycle: maturing
    """

    __slots__ = ()

    def coos(self) -> ReadIter[pa.SparseCOOTensor]:
        raise NotImplementedError

    def dense_tensors(self) -> ReadIter[pa.Tensor]:
        raise NotImplementedError

    def record_batches(self) -> ReadIter[pa.RecordBatch]:
        raise NotImplementedError

    def tables(self) -> ReadIter[pa.Table]:
        raise NotImplementedError
