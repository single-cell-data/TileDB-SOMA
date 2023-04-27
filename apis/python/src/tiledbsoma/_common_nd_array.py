# Copyright (c) 2021-2023 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2023 TileDB, Inc.
#
# Licensed under the MIT License.

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
from ._types import OpenTimestamp
from .options._soma_tiledb_context import (
    SOMATileDBContext,
    _validate_soma_tiledb_context,
)
from .options._tiledb_create_options import TileDBCreateOptions


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
        tiledb_timestamp: Optional[OpenTimestamp] = None,
    ) -> Self:
        """Creates a SOMA ``NDArray`` at the given URI.

        Args:

            type:
                The Arrow type to be stored in the NDArray.
                If the type is unsupported, an error will be raised.
            shape:
                The maximum capacity of each dimension, including room
                for any intended future appends, as a sequence.  E.g. ``(100, 10)``.
                All lengths must be in the postive int64 range, or ``None``.  It's
                necessary to say ``shape=(None, None)`` or ``shape=(None, None,
                None)``, as the sequence length determines the number of dimensions
                N in the N-dimensional array.

                For :class:`SparseNDArray` only, if a slot is None, then the maximum
                possible int64 will be used.  This makes a :class:`SparseNDArray`
                growable.
            platform_config:
                Platform-specific options used to create this Array,
                provided via ``{"tiledb": {"create": ...}}`` nested keys.
            tiledb_timestamp:
                If specified, overrides the default timestamp
                used to open this object. If unset, uses the timestamp provided by
                the context.

        Returns:
            The created NDArray.

        Raises:
            TypeError:
                If the ``type`` is unsupported.
            ValueError:
                If the ``shape`` is unsupported.
            TileDBError:
                If unable to create the underlying object.

        Lifecycle:
            Experimental.
        """
        # Implementor note: we carefully say "maximum possible int64 size" rather than 2**63-1. The
        # reason that the latter, while temptingly simple, is actually untrue is that tiledb core
        # requires that the capacity, when rounded up to an exact multiple of the extent, needs to
        # be representable as a signed 64-bit integer.  So in particular when a unit test (or anyone
        # else) sets extent to a not-power-of-two number like 999 or 1000 then the create fails if
        # the upper limit is exactly 2**63-1.

        context = _validate_soma_tiledb_context(context)
        schema = cls._build_tiledb_schema(
            type,
            shape,
            TileDBCreateOptions.from_platform_config(platform_config),
            context,
            is_sparse=cls.is_sparse,
        )
        handle = cls._create_internal(uri, schema, context, tiledb_timestamp)
        return cls(
            handle,
            _dont_call_this_use_create_or_open_instead="tiledbsoma-internal-code",
        )

    @property
    def shape(self) -> Tuple[int, ...]:
        """Returns capacity of each dimension, always a list of length ``ndim``.

        Lifecycle:
            Experimental.
        """
        return cast(Tuple[int, ...], tuple(self._soma_reader().shape))

    def reshape(self, shape: Tuple[int, ...]) -> None:
        """Unsupported operation for this object type.

        Lifecycle:
            Experimental.
        """
        raise NotImplementedError("reshape operation not implemented.")

    @classmethod
    def _build_tiledb_schema(
        cls,
        type: pa.DataType,
        shape: Sequence[Union[int, None]],
        create_options: TileDBCreateOptions,
        context: SOMATileDBContext,
        *,
        is_sparse: bool,
    ) -> tiledb.ArraySchema:
        _util.check_type("type", type, (pa.DataType,))

        if not pa.types.is_primitive(type):
            raise TypeError(
                f"Unsupported type {type} --"
                " SOMA NDArrays only support primtive Arrow types"
            )

        if not shape:
            raise ValueError("SOMA NDArrays must have a nonzero number of dimensions")

        dims = []
        for dim_idx, dim_shape in enumerate(shape):
            dim_name = f"soma_dim_{dim_idx}"
            dim_capacity, dim_extent = cls._dim_capacity_and_extent(
                dim_name, dim_shape, create_options
            )
            dim = tiledb.Dim(
                name=dim_name,
                domain=(0, dim_capacity - 1),
                tile=dim_extent,
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
            allows_duplicates=create_options.allows_duplicates(),
            offsets_filters=create_options.offsets_filters(),
            validity_filters=create_options.validity_filters(),
            capacity=create_options.capacity(),
            tile_order=tile_order,
            cell_order=cell_order,
            ctx=context.tiledb_ctx,
        )

    @classmethod
    def _dim_capacity_and_extent(
        cls,
        dim_name: str,
        dim_shape: Optional[int],
        create_options: TileDBCreateOptions,
    ) -> Tuple[int, int]:
        raise NotImplementedError("must be implemented by child class.")
