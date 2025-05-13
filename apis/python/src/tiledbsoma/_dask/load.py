# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""Utilities for loading a TileDB-SOMA ``X`` matrix as a Dask Array."""

from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Any,
)

import numpy as np
import numpy.typing as npt
import pyarrow as pa
import scipy.sparse as sp
from somacore.options import PlatformConfig, ResultOrder, ResultOrderStr, SparseNDCoords

from tiledbsoma._dask.util import (
    ChunkSize,
    JoinIDs,
    SOMADaskConfig,
    chunk_ids_sizes,
    coords_to_joinids,
    make_context,
)
from tiledbsoma._dense_nd_array import DenseNDArray
from tiledbsoma._fastercsx import CompressedMatrix, Format
from tiledbsoma._indexer import IntIndexer
from tiledbsoma._sparse_nd_array import SparseNDArray
from tiledbsoma._types import OpenTimestamp
from tiledbsoma.options._soma_tiledb_context import ConfigDict, ConfigVal

if TYPE_CHECKING:
    try:
        import dask.array as da
    except ImportError:
        pass


__all__ = [
    # Re-export `DaskConfig` here, to save most users from needing separate `tiledbsoma._dask` import lines. It also
    # appears in the public API of `SparseNDArrayRead.dask_array`, which imports it directly from
    # `tiledbsoma._dask.util` (to avoid a circular dependency, since `SparseNDArray` is also referenced in this file).
    "SOMADaskConfig",
    "load_daskarray",
    "sparse_chunk",
]


def sparse_chunk(
    block: npt.NDArray[tuple[JoinIDs, JoinIDs]],
    block_info: dict[int | None, Any],
    uri: str,
    tiledb_config: ConfigDict,
    format: Format = "csr",
    result_order: ResultOrderStr = ResultOrder.AUTO,
    tiledb_timestamp: OpenTimestamp | None = None,
    platform_config: PlatformConfig | None = None,
) -> sp.csr_matrix | sp.csc_matrix:
    """Load a slice of a TileDB-SOMA ``X`` matrix, as a block of a CSR- or CSC-backed Dask Array."""
    shape = block_info[None]["chunk-shape"]
    assert block.shape == (1, 1)
    joinids = block[0, 0]
    obs_joinids, var_joinids = joinids
    obs_joinids = obs_joinids.astype("int64")
    var_joinids = var_joinids.astype("int64")
    soma_ctx = make_context(tiledb_config=tiledb_config)
    with SparseNDArray.open(
        uri,
        context=soma_ctx,
        tiledb_timestamp=tiledb_timestamp,
        platform_config=platform_config,
    ) as arr:
        tbls = arr.read(
            coords=(obs_joinids, var_joinids),
            result_order=result_order,
            platform_config=platform_config,
        ).tables()
        # TODO: i32-based reindexing?
        obs_indexer = IntIndexer(obs_joinids, context=soma_ctx)
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
            format=format,
            make_sorted=True,
            context=soma_ctx,
        )
        csx = cm.to_scipy()
        return csx


def load_daskarray(
    layer: SparseNDArray | DenseNDArray,
    *,
    coords: SparseNDCoords | None = None,
    chunk_size: ChunkSize | None = None,
    tiledb_config: dict[str, ConfigVal] | None = None,
    format: Format = "csr",
    result_order: ResultOrderStr = ResultOrder.AUTO,
    platform_config: PlatformConfig | None = None,
) -> "da.Array":
    """Load a TileDB-SOMA X layer as a Dask array."""
    import dask.array as da

    if chunk_size is None:
        raise ValueError("chunk_size required (directly or via `config`)")

    _, _, data_dtype = layer.schema.types
    dtype = data_dtype.to_pandas_dtype()
    nobs, nvar = layer.shape

    obs_ids, var_ids = coords_to_joinids(coords, shape=layer.shape)

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
        if tiledb_config:
            tiledb_config = {**layer.context.tiledb_config}
            tiledb_config.update(**tiledb_config)
        else:
            tiledb_config = layer.context.tiledb_config

        chunk_joinids = da.from_array(arr, chunks=(1, 1))
        meta = (
            sp.csr_matrix((0, 0), dtype=dtype)
            if format == "csr"
            else sp.csc_matrix((0, 0), dtype=dtype)
        )
        X = da.map_blocks(
            sparse_chunk,
            chunk_joinids,
            chunks=(tuple(obs_chunk_sizes), tuple(var_chunk_sizes)),
            meta=meta,
            uri=layer.uri,
            format=format,
            tiledb_config=tiledb_config,
            result_order=result_order,
            tiledb_timestamp=layer.tiledb_timestamp_ms,
            platform_config=platform_config,
        )
    else:
        # TODO: combine with Spatial DenseNDArray Dask code
        raise NotImplementedError(
            f"Dask-loading not implemented yet for DenseNDArray ({layer.uri})"
        )
    return X
