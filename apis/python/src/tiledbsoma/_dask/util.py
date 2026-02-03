# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""Utilities for loading a TileDB-SOMA ``X`` matrix as a Dask Array."""

from __future__ import annotations

from collections.abc import Sequence
from functools import cache
from typing import (
    Any,
    TypedDict,
    Union,
)

import numpy as np
import pyarrow as pa
from numpy import int64
from numpy.typing import NDArray
from typing_extensions import TypeAlias

from tiledbsoma._core_options import SparseNDCoord, SparseNDCoords
from tiledbsoma._soma_context import SOMAContext

ChunkSize = Union[int, tuple[Union[int, None], int]]
JoinIDs: TypeAlias = NDArray[int64]


class SOMADaskConfig(TypedDict, total=False):
    """Dask-related configs.

    Sometimes a ``SOMADaskConfig`` is passed via a ``dask=dict(…)`` kwarg (e.g. to methods that may or may not execute
    in "Dask mode", depending on presence/absence of a ``dask`` config ``dict``). Other methods, that are always
    Dask-related, can receive a "spread" ``SOMADaskConfig`` as ``**kwargs``, ``TypedDict`` supports type-checking for
    both ``kwargs`` styles.
    """

    chunk_size: ChunkSize
    tiledb_config: dict[str, str | float]


def chunk_ids_sizes(joinids: JoinIDs, chunk_size: int, dim_size: int) -> tuple[list[JoinIDs], list[int]]:  # noqa: ARG001
    """Slice chunks from joinids, return chunks' joinids and sizes."""
    chunk_joinids: list[JoinIDs] = []
    chunk_sizes: list[int] = []
    i = 0
    n = len(joinids)
    while i < n:
        end_idx = min(n, i + chunk_size)
        num = end_idx - i
        chunk = joinids[i:end_idx]
        chunk_joinids.append(chunk)
        chunk_sizes.append(num)
        i += num
    return chunk_joinids, chunk_sizes


def make_context(tiledb_config: dict[str, Any]) -> SOMAContext:
    """Create and cache ``SOMATileDBContext``s within Dask worker processes.

    This wrapper just flattens the ``tiledb_config`` ``dict`` into hashable ``tuple``s (for use with ``cache``).
    """
    return _make_context(tiledb_configs=tuple(tiledb_config.items()))


@cache
def _make_context(tiledb_configs: tuple[tuple[str, Any], ...]) -> SOMAContext:
    """Create and cache ``SOMATileDBContext``s within Dask worker processes.

    ``tiledb_config`` is conceptually a ``dict``, but flattened into hashable ``tuple``s here, for use with ``cache``.
    """
    return SOMAContext.create(config=dict(tiledb_configs))


def coord_to_joinids(coord: SparseNDCoord, n: int) -> JoinIDs:
    if not coord:
        return np.arange(n)
    if isinstance(coord, (pa.IntegerArray, pa.ChunkedArray)):
        return coord.to_numpy()
    if isinstance(coord, slice):
        start = coord.start or 0
        stop = coord.stop or n
        step = coord.step or 1
        return np.arange(start, stop, step)
    if isinstance(coord, Sequence):
        return np.array(coord)
    if isinstance(coord, int):
        return np.array([coord])
    raise ValueError(f"Unexpected coord type {type(coord)}: {coord}")


def coords_to_joinids(coords: SparseNDCoords | None, shape: tuple[int, ...]) -> tuple[JoinIDs, JoinIDs]:
    """Convert a ``SparseNDCoords`` to two Numpy arrays (for obs and var), for slicing into Dask tasks."""
    if len(shape) != 2:
        raise ValueError(f"Shape must have length 2: {shape}")
    n_obs, n_var = shape
    if not coords:
        obs = var = None
    elif len(coords) == 1:
        obs = coords[0]
        var = None
    elif len(coords) == 2:
        obs, var = coords
    else:
        raise ValueError(f"coords must be a list of 0, 1, or 2 elements, for {len(coords)}: {coords}")
    return coord_to_joinids(obs, n_obs), coord_to_joinids(var, n_var)
