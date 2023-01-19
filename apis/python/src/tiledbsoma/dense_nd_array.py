from typing import Any, List, Optional, Union, cast

import numpy as np
import pyarrow as pa
import somacore
import tiledb
from somacore import options

# This package's pybind11 code
import tiledbsoma.libtiledbsoma as clib
import tiledbsoma.util as util
import tiledbsoma.util_arrow as util_arrow
from tiledbsoma.util import dense_indices_to_shape

from .collection import CollectionBase
from .exception import SOMAError
from .options import SOMATileDBContext, TileDBCreateOptions
from .tiledb_array import TileDBArray
from .types import NTuple, PlatformConfig

_UNBATCHED = options.BatchSize()


class DenseNDArray(TileDBArray, somacore.DenseNDArray):
    """
    Represents ``X`` and others.
    """

    def __init__(
        self,
        uri: str,
        *,
        parent: Optional[CollectionBase[Any]] = None,
        context: Optional[SOMATileDBContext] = None,
    ):
        """
        Also see the ``TileDBObject`` constructor.
        """
        super().__init__(uri=uri, parent=parent, context=context)

    # Inherited from somacore
    # soma_type: Final = "SOMADenseNDArray"

    def create(
        self,
        type: pa.DataType,
        shape: Union[NTuple, List[int]],
        platform_config: Optional[somacore.options.PlatformConfig] = None,
    ) -> "DenseNDArray":
        """
        Create a ``DenseNDArray`` named with the URI.

        :param type: an Arrow type defining the type of each element in the array. If the type is unsupported, an error will be raised.

        :param shape: the length of each domain as a list, e.g., [100, 10]. All lengths must be in the positive int64 range.

        :param platform_config: Platform-specific options used to create this Array, provided via "tiledb"->"create" nested keys
        """

        # check on shape
        if len(shape) == 0 or any(e <= 0 for e in shape):
            raise ValueError(
                "DenseNDArray shape must be non-zero length tuple of ints > 0"
            )

        if not pa.types.is_primitive(type):
            raise TypeError(
                "Unsupported type - DenseNDArray only supports primtive Arrow types"
            )

        tiledb_create_options = TileDBCreateOptions.from_platform_config(
            platform_config
        )

        dims = []
        for n, e in enumerate(shape):
            dim_name = f"soma_dim_{n}"
            dim = tiledb.Dim(
                name=dim_name,
                domain=(0, e - 1),
                tile=tiledb_create_options.dim_tile(dim_name, min(e, 2048)),
                dtype=np.int64,
                filters=tiledb_create_options.dim_filters(
                    dim_name,
                    [
                        dict(
                            _type="ZstdFilter",
                            level=tiledb_create_options.string_dim_zstd_level(),
                        )
                    ],
                ),
            )
            dims.append(dim)
        dom = tiledb.Domain(dims, ctx=self._ctx)

        attrs = [
            tiledb.Attr(
                name="soma_data",
                dtype=util_arrow.tiledb_type_from_arrow_type(type),
                filters=tiledb_create_options.attr_filters("soma_data", ["ZstdFilter"]),
                ctx=self._ctx,
            )
        ]

        cell_order, tile_order = tiledb_create_options.cell_tile_orders()

        sch = tiledb.ArraySchema(
            domain=dom,
            attrs=attrs,
            sparse=False,
            allows_duplicates=False,
            offsets_filters=tiledb_create_options.offsets_filters(),
            capacity=tiledb_create_options.get("capacity", 100000),
            cell_order=cell_order,
            tile_order=tile_order,
            ctx=self._ctx,
        )

        tiledb.Array.create(self._uri, sch)

        self._common_create(self.soma_type)  # object-type metadata etc

        return self

    @property
    def shape(self) -> NTuple:
        """
        Return length of each dimension, always a list of length ``ndim``
        """
        with self._tiledb_open() as A:
            return cast(NTuple, A.schema.domain.shape)

    def reshape(self, shape: NTuple) -> None:
        """
        Unsupported operation for this object type [lifecycle: experimental].
        """
        raise NotImplementedError("reshape operation not implemented.")

    # Inherited from somacore
    # * ndim accessor
    # * is_sparse: Final = False

    def read(
        self,
        coords: options.DenseNDCoords,
        *,
        result_order: options.StrOr[
            somacore.ResultOrder
        ] = somacore.ResultOrder.ROW_MAJOR,
        batch_size: options.BatchSize = _UNBATCHED,
        partitions: Optional[options.ReadPartitions] = None,
        platform_config: Optional[PlatformConfig] = None,
    ) -> pa.Tensor:
        """
        Read a user-defined dense slice of the array and return as an Arrow ``Tensor``.
        Coordinates must specify a contiguous subarray, and the number of coordinates
        must be less than or equal to the number of dimensions. For example, if the array
        is 10 by 20, then some acceptable values of ``coords`` include ``(3, 4)``,
        ``(slice(5, 10),)``, and ``(slice(5, 10), slice(6, 12))``. Slice indices are
        doubly inclusive.
        """
        del batch_size, partitions, platform_config  # Currently unused.
        result_order = somacore.ResultOrder(result_order)

        with self._tiledb_open("r") as A:
            target_shape = dense_indices_to_shape(coords, A.shape, result_order)
            schema = A.schema
            ned = A.nonempty_domain()

        sr = clib.SOMAReader(
            self._uri,
            name=self.__class__.__name__,
            result_order=result_order.value,
            platform_config=self._ctx.config().dict(),
        )

        if coords is not None:
            if not isinstance(coords, (list, tuple)):
                raise TypeError(
                    f"coords type {type(coords)} unsupported; expected list or tuple"
                )
            if len(coords) < 1 or len(coords) > schema.domain.ndim:
                raise ValueError(
                    f"coords {coords} must have length between 1 and ndim ({schema.domain.ndim}); got {len(coords)}"
                )

            for i, coord in enumerate(coords):
                dim_name = schema.domain.dim(i).name
                if coord is None:
                    pass  # No constraint; select all in this dimension
                elif isinstance(coord, int):
                    sr.set_dim_points(dim_name, [coord])
                elif isinstance(coord, slice):
                    lo_hi = util.slice_to_range(coord, ned[i]) if ned else None
                    if lo_hi is not None:
                        lo, hi = lo_hi
                        if lo < 0 or hi < 0:
                            raise ValueError(
                                f"slice start and stop may not be negative; got ({lo}, {hi})"
                            )
                        if lo > hi:
                            raise ValueError(
                                f"slice start must be <= slice stop; got ({lo}, {hi})"
                            )
                        sr.set_dim_ranges(dim_name, [lo_hi])
                    # Else, no constraint in this slot. This is `slice(None)` which is like
                    # Python indexing syntax `[:]`.
                else:
                    raise TypeError(f"coord type {type(coord)} at slot {i} unsupported")

        sr.submit()

        arrow_tables = []
        while True:
            arrow_table_piece = sr.read_next()
            if not arrow_table_piece:
                break
            arrow_tables.append(arrow_table_piece)

        # For dense arrays there is no zero-output case: attempting to make a test case
        # to do that, say by indexing a 10x20 array by positions 888 and 999, results
        # in read-time errors of the form
        #
        # [TileDB::Subarray] Error: Cannot add range to dimension 'soma_dim_0'; Range [888, 888] is
        # out of domain bounds [0, 9]
        if not arrow_tables:
            raise SOMAError(
                "internal error: at least one table-piece should have been returned"
            )

        arrow_table = pa.concat_tables(arrow_tables)
        return pa.Tensor.from_numpy(
            arrow_table.column("soma_data").to_numpy().reshape(target_shape)
        )

    def write(
        self,
        coords: options.DenseNDCoords,
        values: pa.Tensor,
        *,
        platform_config: Optional[PlatformConfig] = None,
    ) -> None:
        """
        Write subarray, defined by ``coords`` and ``values``. Will overwrite existing
        values in the array.

        Parameters
        ----------
        coords - per-dimension tuple of scalar or slice
            Define the bounds of the subarray to be written.

        values - pyarrow.Tensor
            Define the values to be written to the subarray.  Must have same shape
            as defind by ``coords``, and the type must match the DenseNDArray.
        """
        del platform_config  # Currently unused.
        with self._tiledb_open("w") as A:
            A[coords] = values.to_numpy()
