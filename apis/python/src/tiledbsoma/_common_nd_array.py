"""Common code shared by both NDArray implementations."""

from typing import Optional, Sequence, Tuple, Union, cast

import numpy as np
import pyarrow as pa
import somacore
import tiledb
from somacore import options
from typing_extensions import Self

from . import _arrow_types, _util
from ._tiledb_array import TileDBArray
from .options.soma_tiledb_context import SOMATileDBContext
from .options.tiledb_create_options import TileDBCreateOptions


class NDArray(TileDBArray, somacore.NDArray):
    """Abstract base for the common behaviors of both kinds of NDArray."""

    __slots__ = ()

    @classmethod
    def create(
        cls,
        uri: str,
        *,
        type: pa.DataType,
        shape: Sequence[Union[int, None]],
        platform_config: Optional[options.PlatformConfig] = None,
        context: Optional[SOMATileDBContext] = None,
    ) -> Self:
        """
        Create a SOMA ``NDArray`` at the given URI.

        [lifecycle: experimental]

        :param type: The Arrow type to be stored in the NDArray.
            If the type is unsupported, an error will be raised.

        :param shape: The maximum capacity of the dataframe, including room for any
            intended future appends,as a sequence. E.g. ``(100, 10)``.  All
            lengths must be in the positive int64 range, or `None`.  If a slot
            is `None` -- only supported for `SparseNDArray`, not `DenseNDArray`
            -- then the maximum possible int64 size will be used. This makes a
            `SparseNDArray` growable.

        :param platform_config: Platform-specific options used to create this Array,
            provided via ``{"tiledb": {"create": ...}}`` nested keys.
        """
        # Implementor note: we carefully say "maximum possible in64 size" rather than 2**63-1. The
        # reason that the latter, while temptingly simple, is actually untrue is that tiledb core
        # requires that the capacity, when rounded up to an exact multiple of the extent, needs to
        # be representable as a signed 64-bit integer.  So in particular when a unit test (or anyone
        # else) sets extent to a not-power-of-two number like 999 or 1000 then the create fails if
        # the upper limit is exactly 2**63-1.

        context = context or SOMATileDBContext()
        schema = _build_tiledb_schema(
            type,
            shape,
            TileDBCreateOptions.from_platform_config(platform_config),
            context,
            is_sparse=cls.is_sparse,
        )
        handle = cls._create_internal(uri, schema, context)
        return cls(
            handle,
            _dont_call_this_use_create_or_open_instead="tiledbsoma-internal-code",
        )

    @property
    def shape(self) -> Tuple[int, ...]:
        """
        Return capacity of each dimension, always a list of length ``ndim``
        """
        return cast(Tuple[int, ...], self._handle.schema.domain.shape)

    def reshape(self, shape: Tuple[int, ...]) -> None:
        """
        Unsupported operation for this object type.

        [lifecycle: experimental]
        """
        raise NotImplementedError("reshape operation not implemented.")


def _build_tiledb_schema(
    type: pa.DataType,
    shape: Sequence[Union[int, None]],
    create_options: TileDBCreateOptions,
    context: SOMATileDBContext,
    *,
    is_sparse: bool,
) -> tiledb.ArraySchema:
    _util.check_type("type", type, (pa.DataType,))

    # check on shape
    if is_sparse:
        if len(shape) == 0 or any(e is not None and e <= 0 for e in shape):
            raise ValueError(
                "SOMASparseNDArray shape must be non-zero length tuple of positive ints or Nones"
            )
    else:
        if len(shape) == 0 or any(e is None or e <= 0 for e in shape):
            raise ValueError(
                "SOMADenseNDArray shape must be non-zero length tuple of positive ints"
            )

    if not pa.types.is_primitive(type):
        raise TypeError(
            f"Unsupported type {type} --"
            " SOMA NDArrays only support primtive Arrow types"
        )

    dims = []
    for n, e in enumerate(shape):
        dim_name = f"soma_dim_{n}"

        slot_capacity = 2**63 if e is None else e
        extent = min(slot_capacity, create_options.dim_tile(dim_name, 2048))
        # TileDB requires that each signed-64-bit-int domain slot, rounded up to
        # a multiple of the tile extent in that slot, be representable as a
        # signed 64-bit int. So if the tile extent is 999, say, this will exceed
        # 2**63 - 1.
        if e is None:
            slot_capacity -= extent

        dim = tiledb.Dim(
            name=dim_name,
            domain=(0, slot_capacity - 1),
            tile=extent,
            dtype=np.int64,
            filters=create_options.dim_filters(
                dim_name,
                [
                    dict(
                        _type="ZstdFilter",
                        level=create_options.sparse_nd_array_dim_zstd_level(),
                    )
                ],
            ),
        )
        dims.append(dim)
    dom = tiledb.Domain(dims, ctx=context.tiledb_ctx)

    attrs = [
        tiledb.Attr(
            name="soma_data",
            dtype=_arrow_types.tiledb_type_from_arrow_type(type),
            filters=create_options.attr_filters("soma_data", ["ZstdFilter"]),
            ctx=context.tiledb_ctx,
        )
    ]

    cell_order, tile_order = create_options.cell_tile_orders()

    # TODO: accept more TileDB array-schema options from create_options
    # https://github.com/single-cell-data/TileDB-SOMA/issues/876
    return tiledb.ArraySchema(
        domain=dom,
        attrs=attrs,
        sparse=is_sparse,
        allows_duplicates=create_options.get("allows_duplicates", False),
        offsets_filters=create_options.offsets_filters(),
        validity_filters=create_options.validity_filters(),
        capacity=create_options.get("capacity", 100000),
        tile_order=tile_order,
        cell_order=cell_order,
        ctx=context.tiledb_ctx,
    )
