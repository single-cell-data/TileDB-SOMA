# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

from __future__ import annotations

import datetime
import enum
import pathlib
from collections.abc import Sequence
from typing import TYPE_CHECKING, Any, Callable, TypeVar, Union, get_args

import numpy as np
import numpy.typing as npt
import pandas as pd
import pyarrow as pa
from typing_extensions import Literal, Protocol, TypeGuard

StatusAndReason = tuple[bool, str]
"""Information for whether an upgrade-shape or resize would succeed
if attempted, along with a reason why not."""


if TYPE_CHECKING:
    # `pd.{Series,Index}` require type parameters iff `pandas>=2`. Our pandas dependency (in `setup.py`) is unpinned,
    # which generally resolves to `pandas>=2`, but may be pandas<2 if something else in the user's environment requires
    # that. For type-checking purposes, `.pre-commit-config.yaml` specifies `pandas-stubs>=2`, and we type-check against
    # the `pandas>=2` types here.
    PDSeries = pd.Series[Any]
    PDIndex = pd.Index[Any]

    NPInteger = np.integer[npt.NBitBase]
    NPFloating = np.floating[npt.NBitBase]
    NPNDArray = npt.NDArray[np.number[npt.NBitBase]]
    NPIntArray = npt.NDArray[np.integer[npt.NBitBase]]
    NPIInfo = np.iinfo[NPInteger]
    NPFInfo = np.finfo[NPFloating]
else:
    # When not-type-checking, but running with `pandas>=2`, the "missing" type-params don't affect anything.
    PDSeries = pd.Series
    PDIndex = pd.Index

    # Tests pass `Matrix` (type alias which includes `NPNDArray`, via `DenseMatrix`), as well as other numpy types, to
    # `isinstance`, which causes error "argument 2 cannot be a parameterized generic" using the typedefs in the
    # `TYPE_CHECKING` branch above.
    NPInteger = np.integer
    NPFloating = np.floating
    NPNDArray = np.ndarray
    NPIntArray = np.ndarray
    NPIInfo = np.iinfo
    NPFInfo = np.finfo


Path = Union[str, pathlib.Path]

Ids = Union[list[str], list[bytes], list[int]]

Labels = Union[Sequence[str], PDIndex]

NTuple = tuple[int, ...]

IngestMode = Literal["write", "schema_only", "resume"]  # for static-analysis checks
INGEST_MODES = get_args(IngestMode)  # for run-time checks

# Internal version of ``IngestMode`` that includes "update"; see ``IngestionParams``.
_IngestMode = Union[IngestMode, Literal["update"]]
_INGEST_MODES = (*INGEST_MODES, "update")


OpenTimestamp = Union[int, datetime.datetime]
"""Types that can be used as a timestamp to open a TileDB object.

Integers are treated as milliseconds since the Unix epoch.
"""

ArrowReadResult = Union[
    pa.Table,
    pa.RecordBatch,
    pa.Tensor,
    pa.SparseCOOTensor,
    pa.SparseCSRMatrix,
    pa.SparseCSCMatrix,
]


def is_nonstringy_sequence(it: object) -> bool:
    """Returns true if a sequence is a "normal" sequence and not str or bytes.

    str and bytes are "weird" sequences because iterating them gives you
    another str or bytes instance for each character, and when used as a
    sequence is not what users want.
    """
    return not isinstance(it, (str, bytes)) and isinstance(it, Sequence)


def to_string_tuple(obj: str | Sequence[str]) -> tuple[str, ...]:
    """Returns a tuple of string values.

    If the input is a string, it is returned as a tuple with the string as its
    only item. If it is otherwise a sequence of strings, the sequence is converted
    to a tuple.
    """
    return (obj,) if isinstance(obj, str) else tuple(obj)


def str_or_seq_length(obj: str | Sequence[str]) -> int:
    """Returns the number of str values.

    If input is a string, returns 1. Otherwise, returns the number of strings in the
    sequence.
    """
    return 1 if isinstance(obj, str) else len(obj)


_T = TypeVar("_T")
_T_co = TypeVar("_T_co", covariant=True)


class Slice(Protocol[_T_co]):
    """A slice which stores a certain type of object.

    This protocol describes the built-in ``slice`` type, with a hint to callers
    about what type they should put *inside* the slice.  It is for type
    annotations only and is not runtime-checkable (i.e., you can't do
    ``isinstance(thing, Slice)``), because ``range`` objects also have
    ``start``/``stop``/``step`` and would match, but are *not* slices.
    """

    @property
    def start(self) -> _T_co | None: ...

    @property
    def stop(self) -> _T_co | None: ...

    @property
    def step(self) -> _T_co | None: ...


def is_slice_of(__obj: object, __typ: type[_T]) -> TypeGuard[Slice[_T]]:
    return (
        # We only respect `slice`s proper.
        isinstance(__obj, slice)
        and (__obj.start is None or isinstance(__obj.start, __typ))
        and (__obj.stop is None or isinstance(__obj.stop, __typ))
        and (__obj.step is None or isinstance(__obj.step, __typ))
    )


Metadatum = Union[bytes, float, int, str]
METADATA_TYPES = (bytes, float, int, str)


class SOMABaseTileDBType(enum.Enum):
    SOMAArray = 1
    SOMAGroup = 2


DataProtocol = Literal["tiledbv2", "tiledbv3"]

IntegerArray = Union[npt.NDArray[np.int64], pa.IntegerArray]


class IndexLike(Protocol):
    """The basics of what we expect an Index to be.

    This is a basic description of the parts of the ``pandas.Index`` type
    that we use. It is intended as a rough guide so an implementor can know
    that they are probably passing the right "index" type into a function,
    not as a full specification of the types and behavior of ``get_indexer``.
    """

    def get_indexer(self, target: IntegerArray) -> npt.NDArray[np.intp]:
        """Something compatible with Pandas' Index.get_indexer method."""


IndexFactory = Callable[[IntegerArray], IndexLike]
"""Function that builds an index over the given ``IntegerArray``.

This interface is implemented by the callable ``pandas.Index``.
"""
