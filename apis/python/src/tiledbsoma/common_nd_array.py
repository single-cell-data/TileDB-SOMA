"""Common code shared by both NDArray implementations.

TODO: Extract more stuff.
"""

from typing import Sequence

import numpy as np
import pyarrow as pa
import tiledb

from . import util, util_arrow
from .options.soma_tiledb_context import SOMATileDBContext
from .options.tiledb_create_options import TileDBCreateOptions


def build_tiledb_schema(
    type: pa.DataType,
    shape: Sequence[int],
    create_options: TileDBCreateOptions,
    context: SOMATileDBContext,
    *,
    is_sparse: bool,
) -> tiledb.ArraySchema:
    util.check_type("type", type, (pa.DataType,))

    # check on shape
    if len(shape) == 0 or any(e <= 0 for e in shape):
        raise ValueError("DenseNDArray shape must be non-zero length tuple of ints > 0")

    if not pa.types.is_primitive(type):
        raise TypeError(
            f"Unsupported type {type} --"
            " SOMA NDArrays only support primtive Arrow types"
        )

    dims = []
    for n, e in enumerate(shape):
        dim_name = f"soma_dim_{n}"
        dim = tiledb.Dim(
            name=dim_name,
            domain=(0, e - 1),
            tile=min(e, create_options.dim_tile(dim_name, 2048)),
            dtype=np.int64,
            filters=create_options.dim_filters(
                dim_name,
                [
                    dict(
                        _type="ZstdFilter",
                        level=create_options.string_dim_zstd_level(),
                    )
                ],
            ),
        )
        dims.append(dim)
    dom = tiledb.Domain(dims, ctx=context.tiledb_ctx)

    attrs = [
        tiledb.Attr(
            name="soma_data",
            dtype=util_arrow.tiledb_type_from_arrow_type(type),
            filters=create_options.attr_filters("soma_data", ["ZstdFilter"]),
            ctx=context.tiledb_ctx,
        )
    ]

    cell_order, tile_order = create_options.cell_tile_orders()

    return tiledb.ArraySchema(
        domain=dom,
        attrs=attrs,
        sparse=is_sparse,
        allows_duplicates=is_sparse,
        offsets_filters=create_options.offsets_filters(),
        capacity=create_options.get("capacity", 100000),
        tile_order=tile_order,
        cell_order=cell_order,
        ctx=context.tiledb_ctx,
    )
