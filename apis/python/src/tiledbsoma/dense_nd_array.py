from typing import Optional

import pyarrow as pa
import somacore
from somacore import options
from typing_extensions import Self

from . import util
from .common_nd_array import NDArray
from .exception import SOMAError
from .util import dense_indices_to_shape


class DenseNDArray(NDArray, somacore.DenseNDArray):
    """
    Represents ``X`` and others.

    [lifecycle: experimental]
    """

    __slots__ = ()

    def read(
        self,
        coords: options.DenseNDCoords = (),
        *,
        result_order: options.ResultOrderStr = somacore.ResultOrder.ROW_MAJOR,
        partitions: Optional[options.ReadPartitions] = None,
        platform_config: Optional[options.PlatformConfig] = None,
    ) -> pa.Tensor:
        """
        Read a user-defined dense slice of the array and return as an Arrow ``Tensor``.

        Coordinates must specify a contiguous subarray, and the number of coordinates
        must be less than or equal to the number of dimensions. For example,
        if the array is 10 by 20, then some acceptable values of ``coords`` include
        ``(3, 4)``, ``(slice(5, 10),)``, and ``(slice(5, 10), slice(6, 12))``.
        Slice indices are doubly inclusive.

        [lifecycle: experimental]
        """
        del partitions, platform_config  # Currently unused.
        self._check_open_read()
        result_order = somacore.ResultOrder(result_order)

        schema = self._handle.schema
        target_shape = dense_indices_to_shape(coords, schema.shape, result_order)

        sr = self._soma_reader(result_order=result_order.value)

        if not isinstance(coords, (list, tuple)):
            raise TypeError(
                f"coords type {type(coords)} unsupported; expected list or tuple"
            )
        if len(coords) > schema.domain.ndim:
            raise ValueError(
                f"coords {coords} must be shorter than ndim ({schema.domain.ndim}); got {len(coords)}"
            )

        for i, coord in enumerate(coords):
            dim_name = schema.domain.dim(i).name
            if coord is None:
                pass  # No constraint; select all in this dimension
            elif isinstance(coord, int):
                sr.set_dim_points(dim_name, [coord])
            elif isinstance(coord, slice):
                ned = self._handle.reader.nonempty_domain()
                # ned is None iff the array has no data
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
        platform_config: Optional[options.PlatformConfig] = None,
    ) -> Self:
        """
        Write subarray, defined by ``coords`` and ``values``. Will overwrite existing
        values in the array.

        [lifecycle: experimental]

        :param coords: A per-dimension tuple of scalars or slices
            defining the bounds of the subarray to be written.
        :param values: The values to be written to the subarray.  Must have
        the same shape as ``coords``, and the type must match the DenseNDArray.
        :param platform_config: Optional platform-specific options to use
            in this write operation (currently unused).
        """
        util.check_type("values", values, (pa.Tensor,))

        del platform_config  # Currently unused.
        self._handle.writer[coords] = values.to_numpy()
        return self
