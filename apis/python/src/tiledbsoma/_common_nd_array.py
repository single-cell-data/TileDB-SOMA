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
                possible int32 will be used.  This makes a :class:`SparseNDArray`
                growable.
            platform_config:
                Platform-specific options used to create this array.
                This may be provided as settings in a dictionary, with options
                located in the ``{'tiledb': {'create': ...}}`` key,
                or as a :class:`~tiledbsoma.TileDBCreateOptions` object.
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
        This will not necessarily match the bounds of occupied cells within the array.
        Rather, it is the bounds outside of which no data may be written.

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
                filters=create_options.dim_filters_tiledb(
                    dim_name,
                    [
                        dict(
                            _type="ZstdFilter",
                            level=create_options.sparse_nd_array_dim_zstd_level,
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
                filters=create_options.attr_filters_tiledb("soma_data", ["ZstdFilter"]),
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
            allows_duplicates=create_options.allows_duplicates,
            offsets_filters=create_options.offsets_filters_tiledb(),
            validity_filters=create_options.validity_filters_tiledb(),
            capacity=create_options.capacity,
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
