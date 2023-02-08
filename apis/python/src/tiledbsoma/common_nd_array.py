"""Common code shared by both NDArray implementations."""

from typing import Optional, Sequence, Tuple, cast

import numpy as np
import pyarrow as pa
import somacore
import tiledb
from somacore import options
from typing_extensions import Self

from . import arrow_types, util
from .options.soma_tiledb_context import SOMATileDBContext
from .options.tiledb_create_options import TileDBCreateOptions
from .tiledb_array import TileDBArray


class NDArray(TileDBArray, somacore.NDArray):
    """Abstract base for the common behaviors of both kinds of NDArray."""

    __slots__ = ()

    @classmethod
    def create(
        cls,
        uri: str,
        *,
        type: pa.DataType,
        shape: Sequence[int],
        platform_config: Optional[options.PlatformConfig] = None,
        context: Optional[SOMATileDBContext] = None,
    ) -> Self:
        """
        Create a SOMA ``NDArray`` at the given URI.

        [lifecycle: experimental]

        :param type: The Arrow type to be stored in the NDArray.
            If the type is unsupported, an error will be raised.

        :param shape: the length of each dimension as a sequence, e.g., ``[100, 10]``.
            All lengths must be in the positive int64 range.

        :param platform_config: Platform-specific options used to create this Array,
            provided via ``{"tiledb": {"create": ...}}`` nested keys.
        """
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
        Return length of each dimension, always a list of length ``ndim``
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
            dtype=arrow_types.tiledb_type_from_arrow_type(type),
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
