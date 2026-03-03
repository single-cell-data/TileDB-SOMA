# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.
"""Enums and other types used as options across methods of many types.

These types are *concrete* and should be used as-is as inputs to the various
SOMA types that require them, not reimplemented by the implementing package.
"""

from __future__ import annotations

import enum
from collections.abc import Mapping, Sequence
from typing import Any, Final, TypeVar, Union

import attrs
import numpy as np
import numpy.typing as npt
import pyarrow as pa
import shapely
from typing_extensions import Literal

from ._types import Slice

SOMA_JOINID: Final = "soma_joinid"
"""Global constant for the SOMA join ID."""

SOMA_GEOMETRY: Final = "soma_geometry"
"""Global constant for SOMA spatial geometry type."""

OpenMode = Literal["r", "w", "d"]
"""How to open a SOMA object: read, write, or delete."""


class ReadPartitions:
    """Sentinel base class for read-partition types.

    Lifecycle:
        Experimental
    """

    __slots__ = ()


@attrs.define(frozen=True)
class IOfN(ReadPartitions):
    """Specifies that a read should return partition ``i`` out of ``n`` total.

    For a read operation that returns ``n`` partitions, the read operation will
    return the ``i``th partition (zero-indexed) out of ``n`` partitions of
    approximately equal size.

    Lifecycle:
        Experimental
    """

    i: int = attrs.field()
    """Which partition to return (zero-indexed)."""
    n: int = attrs.field()
    """How many partitions there will be."""

    @i.validator
    def _validate(self, _, __):  # type: ignore[no-untyped-def]  # noqa: ANN202, ANN001
        del _, __  # Unused.
        if not 0 <= self.i < self.n:
            raise ValueError(f"Partition index {self.i} must be in the range [0, {self.n})")


@attrs.define(frozen=True)
class BatchSize:
    """Specifies the size of a batch that should be returned from reads.

    Read operations on foundational types return an iterator over "batches" of
    data, enabling processing of larger-than-core datasets. This class allows
    you to control what the size of those batches is.

    If none of these options are set, a "reasonable" batch size is determined
    automatically.

    For example::

        BatchSize(count=100)
        # Will return batches of 100 elements.

        BatchSize(bytes=1024 ** 2)
        # Will return batches of up to 1 MB.

        BatchSize()
        # Will return automatically-sized batches.

    Lifecycle:
        Experimental
    """

    count: int | None = attrs.field(default=None)
    """``arrow.Table``s with this number of rows will be returned."""
    bytes: int | None = attrs.field(default=None)
    """Data of up to this size in bytes will be returned."""

    @count.validator
    @bytes.validator
    def _validate(self, attr, value: int | None) -> None:  # type: ignore[no-untyped-def]  # noqa: ANN001
        if not value:
            return  # None (or 0, which we treat equivalently) is always valid.
        if value < 0:
            raise ValueError(f"If set, '{attr.name}' must be positive")
        if self.count and self.bytes:
            raise ValueError("Either 'count' or 'bytes' may be set, not both")


PlatformConfig = Union[dict[str, Mapping[str, Any]], object]
"""Type alias for the ``platform_config`` parameter.

``platform_config`` allows platform-specific configuration data to be passed
to individual calls. This is either a ``dict``, or an implementation-specific
configuration object:

- If a dictionary, the keys to the dictionary are the name of a SOMA
  implementation, and the value of the dictionary is configuration specific to
  that implementation.
- If an implementation-specific object, that implementation will use that object
  for configuration data; it will be ignored by others.

See the "Per-call configuration" section of the main SOMA specification.
"""


class ResultOrder(enum.Enum):
    """The order results should be returned in.

    Lifecycle:
        Experimental
    """

    AUTO = "auto"
    ROW_MAJOR = "row-major"
    COLUMN_MAJOR = "column-major"


ResultOrderStr = Union[ResultOrder, Literal["auto", "row-major", "column-major"]]
"""A ResultOrder, or the str representing it."""

DenseCoord = Union[int, Slice[int], None]
"""A single coordinate range for reading dense data.

``None`` indicates the entire domain of a dimension; values of this type are
not ``Optional``, but may be ``None``.
"""


DenseNDCoords = Sequence[DenseCoord]
"""A sequence of ranges to read dense data."""

_T = TypeVar("_T")
ValSliceOrSequence = Union[_T, Slice[_T], Sequence[_T]]
"""A value of a type, a Slice of that type, or a Sequence of that type."""

# NOTE: Keep this in sync with the types accepted in `_canonicalize_coord`
# in ./query/axis.py.
SparseDFCoord = Union[
    ValSliceOrSequence[bytes],
    ValSliceOrSequence[float],
    ValSliceOrSequence[int],
    ValSliceOrSequence[slice],
    ValSliceOrSequence[str],
    ValSliceOrSequence[np.datetime64],
    ValSliceOrSequence[pa.TimestampType],
    pa.Array,
    pa.ChunkedArray,
    npt.NDArray[np.integer],
    npt.NDArray[np.datetime64],
    None,
]
"""A single coordinate range for one dimension of a sparse dataframe."""
SparseDFCoords = Sequence[SparseDFCoord]
"""A sequence of coordinate ranges for reading dense dataframes."""

SparseNDCoord = Union[
    ValSliceOrSequence[int],
    npt.NDArray[np.integer],
    pa.IntegerArray,
    pa.ChunkedArray,
    None,
]
"""A single coordinate range for one dimension of a sparse ndarray."""

SparseNDCoords = Sequence[SparseNDCoord]
"""A sequence of coordinate ranges for reading sparse ndarrays."""

SpatialRegion = Union[Sequence[int], Sequence[float], shapely.geometry.base.BaseGeometry]
"""A spatial region used for reading spatial dataframes and multiscale images."""
