from typing import Any, List, Literal, Optional, Tuple, Union, cast

import numpy as np
import pyarrow as pa
import tiledb

import tiledbsoma.util_arrow as util_arrow
from tiledbsoma.util_tiledb import tiledb_result_order_from_soma_result_order

from .soma_collection import CollectionBase
from .tiledb_array import TileDBArray
from .tiledb_platform_config import TileDBPlatformConfig
from .types import NTuple, SOMADenseNdCoordinates, SOMAResultOrder


class SOMADenseNdArray(TileDBArray):
    """
    Represents ``X`` and others.
    """

    def __init__(
        self,
        uri: str,
        *,
        parent: Optional[CollectionBase[Any]] = None,
        tiledb_platform_config: Optional[TileDBPlatformConfig] = None,
        ctx: Optional[tiledb.Ctx] = None,
    ):
        """
        Also see the ``TileDBObject`` constructor.
        """
        super().__init__(
            uri=uri,
            parent=parent,
            tiledb_platform_config=tiledb_platform_config,
            ctx=ctx,
        )

    @property
    def soma_type(self) -> Literal["SOMADenseNdArray"]:
        return "SOMADenseNdArray"

    def create(
        self,
        type: pa.DataType,
        shape: Union[NTuple, List[int]],
    ) -> "SOMADenseNdArray":
        """
        Create a ``SOMADenseNdArray`` named with the URI.

        :param type: an Arrow type defining the type of each element in the array. If the type is unsupported, an error will be raised.

        :param shape: the length of each domain as a list, e.g., [100, 10]. All lengths must be in the positive int64 range.
        """

        # check on shape
        if len(shape) == 0 or any(e <= 0 for e in shape):
            raise ValueError(
                "DenseNdArray shape must be non-zero length tuple of ints > 0"
            )

        if not pa.types.is_primitive(type):
            raise TypeError(
                "Unsupported type - SOMADenseNdArray only supports primtive Arrow types"
            )

        level = self._tiledb_platform_config.string_dim_zstd_level

        dims = []
        for n, e in enumerate(shape):
            dim = tiledb.Dim(
                name=f"soma_dim_{n}",
                domain=(0, e - 1),
                tile=min(e, 2048),  # TODO: PARAMETERIZE
                dtype=np.int64,
                filters=[tiledb.ZstdFilter(level=level)],
            )
            dims.append(dim)
        dom = tiledb.Domain(dims, ctx=self._ctx)

        attrs = [
            tiledb.Attr(
                name="soma_data",
                dtype=util_arrow.tiledb_type_from_arrow_type(type),
                filters=[tiledb.ZstdFilter()],
                ctx=self._ctx,
            )
        ]

        sch = tiledb.ArraySchema(
            domain=dom,
            attrs=attrs,
            sparse=False,
            allows_duplicates=False,
            offsets_filters=[
                tiledb.DoubleDeltaFilter(),
                tiledb.BitWidthReductionFilter(),
                tiledb.ZstdFilter(),
            ],
            capacity=100000,
            cell_order="row-major",
            tile_order="row-major",
            ctx=self._ctx,
        )

        tiledb.Array.create(self._uri, sch, ctx=self._ctx)

        self._common_create()  # object-type metadata etc

        return self

    @property
    def shape(self) -> NTuple:
        """
        Return length of each dimension, always a list of length ``ndim``
        """
        with self._tiledb_open() as A:
            return cast(NTuple, A.schema.domain.shape)

    def reshape(self, shape: NTuple) -> None:
        raise NotImplementedError("reshape operation not implemented.")

    @property
    def ndim(self) -> int:
        """
        Return number of index columns
        """
        with self._tiledb_open() as A:
            return cast(int, A.schema.domain.ndim)

    @property
    def is_sparse(self) -> Literal[False]:
        """
        Returns ``False``.
        """
        return False

    def read_tensor(
        self,
        coords: SOMADenseNdCoordinates,
        *,
        result_order: SOMAResultOrder = "row-major",
    ) -> pa.Tensor:
        """
        Read a user-defined dense slice of the array and return as an Arrow ``Tensor``.
        """
        tiledb_result_order = tiledb_result_order_from_soma_result_order(
            result_order, accept=["column-major", "row-major"]
        )
        with self._tiledb_open("r") as A:
            target_shape = _dense_index_to_shape(coords, A.shape, result_order)
            query = A.query(return_arrow=True, order=tiledb_result_order)
            arrow_tbl = query.df[coords]
            return pa.Tensor.from_numpy(
                arrow_tbl.column("soma_data").to_numpy().reshape(target_shape)
            )

    def read_numpy(
        self,
        coords: SOMADenseNdCoordinates,
        *,
        result_order: SOMAResultOrder = "row-major",
    ) -> np.ndarray:
        """
        Read a user-specified dense slice of the array and return as an Numpy ``ndarray``.
        """
        return cast(
            np.ndarray, self.read_tensor(coords, result_order=result_order).to_numpy()
        )

    def write_tensor(
        self,
        coords: SOMADenseNdCoordinates,
        values: pa.Tensor,
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
            as defind by ``coords``, and the type must match the SOMADenseNdArray.
        """
        with self._tiledb_open("w") as A:
            A[coords] = values.to_numpy()

    def write_numpy(self, coords: SOMADenseNdCoordinates, values: np.ndarray) -> None:
        """ "
        Write a numpy ``ndarray`` to the user specified coordinates
        """
        self.write_tensor(coords, pa.Tensor.from_numpy(values))


# module-private utility
def _dense_index_to_shape(
    coords: Tuple[Union[int, slice], ...],
    array_shape: Tuple[int, ...],
    result_order: SOMAResultOrder,
) -> Tuple[int, ...]:
    """
    Given a subarray index specified as a tuple of per-dimension slices or scalars
    (eg, ``([:], 1, [1:2])``), and the shape of the array, return the shape of
    the subarray.

    See read_tensor for usage.
    """
    shape: List[int] = []
    for n, idx in enumerate(coords):
        if type(idx) is int:
            shape.append(1)
        elif type(idx) is slice:
            start, stop, step = idx.indices(array_shape[n])
            if step != 1:
                raise ValueError("stepped slice ranges are not supported")
            shape.append(stop - start)
        else:
            raise ValueError("coordinates must be tuple of int or slice")

    if result_order == "row-major":
        return tuple(shape)

    return tuple(reversed(shape))
