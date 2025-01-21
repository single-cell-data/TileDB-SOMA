# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

from __future__ import annotations

import json
import warnings
from typing import Any, Mapping, Sequence, Tuple, Union

import dask.array as da
import numpy as np

from ._util import _version_less_than

try:
    import spatialdata as sd
    from spatialdata.models.models import DataTree
except ImportError as err:
    warnings.warn("Experimental spatial exporter requires the spatialdata package.")
    raise err
try:
    import xarray as xr
except ImportError as err:
    warnings.warn("Experimental spatial exporter requires the xarray package.")
    raise err

from ... import DenseNDArray
from ...options._soma_tiledb_context import SOMATileDBContext
from ._util import _str_to_int


class DenseNDArrayWrapper:
    """Wraps an SOMA DenseNDArray for array-like accessors from other libraries.


    Args:
        uri: Location of the :class:`DenseNDArray`.
        context: Optional :class:`SOMATileDBContext` containing storage parameters, etc.
    """

    def __init__(
        self,
        uri: str,
        *,
        context: SOMATileDBContext | None = None,
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
        return result.reshape(output_shape)

    @property
    def ndim(self) -> int:
        """Returns the number of dimensions."""
        return len(self._array.shape)

    @property
    def dtype(self) -> np.typing.DTypeLike:
        """Numpy dtype of the array data."""
        return self._dtype

    def recommend_chunks(self) -> Tuple[int, ...]:
        """Returns recommended chunk sizes for chunking this array."""
        dim_info = json.loads(self._array.schema_config_options().dims)
        return tuple(
            _str_to_int(dim_info[f"soma_dim_{index}"]["tile"])
            for index in range(self.ndim)
        )

    @property
    def shape(self) -> Tuple[int, ...]:
        """Shape of the wrapped SOMA DenseNDArray."""
        return self._array.shape


def dense_nd_array_to_data_array(
    uri: str,
    *,
    dim_names: Tuple[str, ...],
    chunks: Tuple[int, ...] | None = None,
    attrs: Mapping[str, Any] | None = None,
    context: SOMATileDBContext | None = None,
) -> xr.DataArray:
    """Create a :class:`xarray.DataArray` that accesses a SOMA :class:`DenseNDarray`
    through dask.

    Args:
        uri: The location of the array to open,.
        dim_names: The name of the dimensions to use for the :class:`xarray.DataArray`.
        chunks: If provided, the chunk sizes to use for the dask array wrapping the
            SOMA :class:`DenseNDArray`. If not provided, the tile size of the array
            will be used.
        attrs: The attributes (metadata) to set on the :class:`xarray.DataArray`.
        context: Optional :class:`SOMATileDBContext` containing storage parameters, etc.
    """

    array_wrapper = DenseNDArrayWrapper(uri=uri, context=context)

    if chunks is None:
        chunks = array_wrapper.recommend_chunks()

    data = da.from_array(
        array_wrapper,
        name="data",
        chunks=chunks,
        asarray=True,
        fancy=False,
    )

    return xr.DataArray(data, dims=dim_names, attrs=attrs)


def images_to_datatree(image_data_arrays: Sequence[xr.DataArray]) -> DataTree:
    # If SpatialData version < 0.2.6 use the legacy xarray_datatree implementation
    # of the DataTree.
    if _version_less_than(sd.__version__, (0, 2, 5)):
        return DataTree.from_dict(
            {f"scale{index}": image for index, image in enumerate(image_data_arrays)}
        )
    return DataTree.from_dict(
        {
            f"scale{index}": xr.Dataset({"image": image})
            for index, image in enumerate(image_data_arrays)
        }
    )
