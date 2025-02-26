# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""Utilities for loading a TileDB-SOMA ``X`` matrix as a Dask Array.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from typing import (
    TYPE_CHECKING,
    Any,
    Union,
)

import numpy as np
import numpy.typing as npt
import pyarrow as pa
import scipy.sparse as sp
from somacore import DenseNDArray

from tiledbsoma._indexer import IntIndexer

from ._sparse_nd_array import SparseNDArray
from .options import SOMATileDBContext
from .options._soma_tiledb_context import ConfigVal

if TYPE_CHECKING:
    from ._experiment import Experiment  # noqa


from functools import cache

if TYPE_CHECKING:
    pass
from ._fastercsx import CompressedMatrix

if TYPE_CHECKING:
    try:
        import dask.array as da
    except ImportError:
        pass


ChunkSize = Union[int, tuple[Union[int, None], int]]


@dataclass
class DaskConfig:
    chunk_size: ChunkSize
    tdb_concurrency: int | None = 1
    tdb_configs: dict[str, ConfigVal] = field(default_factory=dict)


def chunk_ids_sizes(
    joinids: list[int], chunk_size: int, dim_size: int
) -> tuple[list[list[int]], list[int]]:
    """Slice chunks from joinids, return chunks' joinids and sizes."""
    chunk_joinids: list[list[int]] = []
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


def make_context(
    tdb_concurrency: int | None,
    tiledb_config: dict[str, Any],
) -> SOMATileDBContext:
    """Create and cache a ``SOMATileDBContext``, optionally setting several TileDB concurrency configs.

    Intended for use with ``load_daskarray``, where we want to restrict TileDB concurrency when running TileDB-SOMA
    queries from many Dask workers in parallel. In such cases we don't want each TileDB-SOMA query thinking it should
    utilize all ``cpu_count`` CPUs.

    This wrapper just flattens the ``tiledb_config`` ``dict`` into hashable ``tuple``s (for use with ``cache``).
    """
    return _make_context(
        tdb_concurrency=tdb_concurrency,
        tiledb_configs=tuple(tiledb_config.items()),
    )


@cache
def _make_context(
    tdb_concurrency: int | None,
    tiledb_configs: tuple[tuple[str, Any], ...],
) -> SOMATileDBContext:
    """Create and cache a ``SOMATileDBContext``, optionally setting several TileDB concurrency configs.

    Intended for use with ``load_daskarray``, where we want to restrict TileDB concurrency when running TileDB-SOMA
    queries from many Dask workers in parallel. In such cases we don't want each TileDB-SOMA query thinking it should
    utilize all ``cpu_count`` CPUs.

    ``tiledb_config`` is conceptually a ``dict``, but flattened into hashable ``tuple``s here, for use with ``cache``.
    """
    if tdb_concurrency is None:
        tdb_concurrency = 1

    tiledb_config = dict(tiledb_configs)
    threadpool = None
    if tdb_concurrency:
        # `0` omits these kwargs altogether
        tiledb_config.update(
            {
                "sm.io_concurrency_level": tdb_concurrency,
                "sm.compute_concurrency_level": tdb_concurrency,
            }
        )
        threadpool = ThreadPoolExecutor(max_workers=tdb_concurrency)

    return SOMATileDBContext(
        tiledb_config=tiledb_config,
        threadpool=threadpool,
    )


def sparse_chunk(
    block: npt.NDArray[tuple[list[int], list[int]]],
    block_info: dict[int | None, Any],
    uri: str,
    tiledb_config: dict[str, Union[str, float]],
    tdb_concurrency: int | None = None,
) -> sp.csr_matrix:
    """Load a slice of a TileDB-SOMA ``X`` matrix, as a block of a CSR-backed Dask Array."""
    shape = block_info[None]["chunk-shape"]
    assert block.shape == (1, 1)
    joinids = block[0, 0]
    obs_joinids, var_joinids = joinids
    # TODO: create sequence of tables for CM.from_soma; use record batches (idiomatic Arrow)

    soma_ctx = make_context(
        tiledb_config=tiledb_config,
        tdb_concurrency=tdb_concurrency,
    )
    with SparseNDArray.open(uri, context=soma_ctx) as arr:
        tbls = arr.read((obs_joinids, var_joinids)).tables()
        obs_indexer = IntIndexer(
            obs_joinids, context=soma_ctx
        )  # TODO: i32-based reindexing?
        var_indexer = IntIndexer(var_joinids, context=soma_ctx)
        new_tbls = []
        for tbl in tbls:
            new_dim0 = obs_indexer.get_indexer(tbl["soma_dim_0"]).astype("int32")
            new_dim1 = var_indexer.get_indexer(tbl["soma_dim_1"]).astype("int32")
            new_tbl = pa.Table.from_pydict(
                {
                    "soma_dim_0": new_dim0,
                    "soma_dim_1": new_dim1,
                    "soma_data": tbl["soma_data"],
                }
            )
            new_tbls.append(new_tbl)
        cm = CompressedMatrix.from_soma(
            new_tbls,
            shape=shape,
            format="csr",
            make_sorted=True,
            context=soma_ctx,
        )
        csr = cm.to_scipy()
        return csr


def load_daskarray(
    layer: Union[SparseNDArray, DenseNDArray],
    chunk_size: ChunkSize,
    obs_joinids: pa.IntegerArray | None = None,
    var_joinids: pa.IntegerArray | None = None,
    tdb_concurrency: int | None = None,
    **tdb_configs: ConfigVal,
) -> "da.Array":
    """Load a TileDB-SOMA X layer as a Dask array."""
    import dask.array as da

    _, _, data_dtype = layer.schema.types
    dtype = data_dtype.to_pandas_dtype()
    nobs, nvar = layer.shape

    obs_ids = obs_joinids.to_numpy().tolist() if obs_joinids else list(range(nobs))
    var_ids = var_joinids.to_numpy().tolist() if var_joinids else list(range(nvar))

    if isinstance(chunk_size, int):
        obs_chunk_size = chunk_size
        var_chunk_size = nvar
    else:
        _obs_chunk_size, var_chunk_size = chunk_size
        obs_chunk_size = _obs_chunk_size or nobs

    obs_chunk_joinids, obs_chunk_sizes = chunk_ids_sizes(obs_ids, obs_chunk_size, nobs)
    var_chunk_joinids, var_chunk_sizes = chunk_ids_sizes(var_ids, var_chunk_size, nvar)

    arr = np.empty((len(obs_chunk_joinids), len(var_chunk_joinids)), dtype=object)
    for obs_chunk_idx, obs_chunk_ids in enumerate(obs_chunk_joinids):
        for var_chunk_idx, var_chunk_ids in enumerate(var_chunk_joinids):
            arr[obs_chunk_idx, var_chunk_idx] = (obs_chunk_ids, var_chunk_ids)

    if isinstance(layer, SparseNDArray):
        tiledb_config = {**layer.context.tiledb_config}
        tiledb_config.update(**tdb_configs)
        chunk_joinids = da.from_array(arr, chunks=(1, 1))
        X = da.map_blocks(
            sparse_chunk,
            chunk_joinids,
            chunks=(tuple(obs_chunk_sizes), tuple(var_chunk_sizes)),
            meta=sp.csr_matrix((0, 0), dtype=dtype),
            uri=layer.uri,
            tdb_concurrency=tdb_concurrency,
            tiledb_config=tiledb_config,
        )
    else:
        raise NotImplementedError(
            f"Dask-loading not implemented yet for DenseNDArray ({layer.uri})"
        )
    return X
