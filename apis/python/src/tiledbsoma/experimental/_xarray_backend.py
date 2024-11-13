# Copyright (c) 2024 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2024 TileDB, Inc
#
# Licensed under the MIT License.
from typing import Any, Mapping, Optional, Tuple, Union

import dask.array as da
import numpy as np
from xarray import DataArray

from .. import DenseNDArray
from ..options._soma_tiledb_context import SOMATileDBContext


class DenseNDArrayWrapper:
    """Wraps an SOMA DenseNDArray for array-like accessors from other libraries.


    Args:
        uri: Location of the :class:`DenseNDArray`.
    """

    def __init__(
        self,
        uri: str,
        *,
        context: Optional[SOMATileDBContext] = None,
    ):
        self._array = DenseNDArray.open(uri, context=context)
        self._dtype: np.typing.DTypeLike = self._array.schema.field(
            "soma_data"
        ).type.to_pandas_dtype()

    def __getitem__(self, key: Tuple[Union[slice, int], ...]) -> np.typing.NDArray[Any]:
        """Returns a numpy array containing data from a requested tuple."""

        # Compute the expected Xarray output shape.
        output_shape = tuple(
            len(range(dim_size)[index])  # type: ignore
            for dim_size, index in zip(self.shape, key)
            if not np.isscalar(index)
        )

        def update_int_index(index: int, dim_size: int) -> int:
            """Convert xarray integer index to a SOMA-compatible integer index.

            - Convert negative indices to appropriate positive position.
            """
            if not -dim_size <= index < dim_size:
                raise IndexError(
                    f"Index {index} out of bounds for dimension with size {dim_size}."
                )
            return index if index >= 0 else index + dim_size - 1

        def update_slice_index(index: slice, dim_size: int) -> slice:
            """Convert xarray slice index to a SOMA-compatible slice index.

            - Throw error for slice step != 1.
            - Convert negative indices to appropriate positive position.
            - Convert upper index to be inclusive (SOMA/TileDB) instead of
              exclusive (numpy/xarray).
            """
            if index.step not in (1, None):
                raise NotImplementedError("Slice steps are not supported.")
            _index = range(dim_size)[index]  # Convert negative values to positive.
            return slice(_index.start, _index.stop - 1)

        key = tuple(
            (
                update_slice_index(index, dim_size)
                if isinstance(index, slice)
                else update_int_index(index, dim_size)
            )
            for index, dim_size in zip(key, self.shape)
        )

        # Read the data from SOMA, convert to numpy, and reshape.
        result = self._array.read(key).to_numpy()
        return result.reshape(output_shape)  # type: ignore

    @property
    def chunks(self) -> Tuple[int, ...]:
        """Recommended chunk sizes for chunking this array."""
        # TODO: Fix this to return tile sizes instead.
        return self.shape

    @property
    def ndim(self) -> int:
        """Returns the number of dimensions."""
        return len(self._array.shape)

    @property
    def dtype(self) -> np.typing.DTypeLike:
        """Numpy dtype of the array data."""
        return self._dtype

    @property
    def shape(self) -> Tuple[int, ...]:
        """Shape of the wrapped SOMA DenseNDArray."""
        return self._array.shape


def dense_nd_array_to_data_array(
    uri: str,
    *,
    dim_names: Tuple[str, ...],
    chunks: Optional[Tuple[int, ...]] = None,
    attrs: Optional[Mapping[str, Any]] = None,
    context: Optional[SOMATileDBContext] = None,
) -> DataArray:
    """Create a :class:`xarray.DataArray` that accesses a SOMA :class:`DenseNDarray`
    through dask.

    Args:

    """

    array_wrapper = DenseNDArrayWrapper(uri=uri, context=context)

    if chunks is None:
        chunks = array_wrapper.chunks

    data = da.from_array(
        array_wrapper,
        name="data",
        chunks=chunks,
        asarray=True,
        fancy=False,
    )

    return DataArray(data, dims=dim_names, attrs=attrs)
