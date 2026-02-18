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
from dataclasses import dataclass
from typing import Any, Generic, TypeVar, Union

import pyarrow as pa

from ._coordinate_space import CoordinateSpace, CoordinateTransform
from ._core_options import BatchSize, ResultOrder

_RO_AUTO = ResultOrder.AUTO
_ReadData = TypeVar("_ReadData")

AxisDomain = Union[tuple[Any, Any], list[Any], None]
Domain = Sequence[AxisDomain]

_UNBATCHED = BatchSize()

SparseArrowData = Union[
    pa.SparseCSCMatrix,
    pa.SparseCSRMatrix,
    pa.SparseCOOTensor,
    pa.Table,
]
"""Any of the sparse data storages provided by Arrow."""


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

        Lifecycle:
            Maturing.
        """
        raise NotImplementedError


class SparseRead:
    """Intermediate type to choose result format when reading a sparse array.

    A query may not be able to return all of these formats. The concrete result
    may raise a ``NotImplementedError`` or may choose to raise a different
    exception (likely a ``TypeError``) containing more specific information
    about why the given format is not supported.

    Lifecycle:
        Maturing.
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


@dataclass
class SpatialRead(Generic[_ReadData]):
    """Reader for spatial data.

    Args:
        data: The data accessor.
        data_coordinate_space: The coordinate space the read data is defined on.
        output_coordinate_space: The requested output coordinate space.
        coordinate_transform: A coordinate transform from the data coordinate space to
            the desired output coordinate space.

    Lifecycle:
        Experimental.
    """

    data: _ReadData
    data_coordinate_space: CoordinateSpace
    output_coordinate_space: CoordinateSpace
    coordinate_transform: CoordinateTransform

    def __post_init__(self) -> None:
        if self.data_coordinate_space.axis_names != self.coordinate_transform.input_axes:
            raise ValueError("Input coordinate transform axis names do not match the data coordinate space.")
        if self.output_coordinate_space.axis_names != self.coordinate_transform.output_axes:
            raise ValueError("Output coordinate transform axis names do not match the output coordinate space.")
