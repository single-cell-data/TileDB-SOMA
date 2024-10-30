# Copyright (c) 2024 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2024 TileDB, Inc
#
# Licensed under the MIT License.
from typing import Any, Mapping, Optional, Tuple, Union

import numpy as np
from somacore import options
from xarray import Variable
from xarray.backends import AbstractDataStore, BackendArray
from xarray.core.indexing import (
    BasicIndexer,
    ExplicitIndexer,
    IndexingSupport,
    LazilyIndexedArray,
    OuterIndexer,
    VectorizedIndexer,
    explicit_indexing_adapter,
)
from xarray.core.utils import Frozen

from .. import DenseNDArray
from ..options._soma_tiledb_context import SOMATileDBContext


class DenseNDArrayWrapper(BackendArray):  # type: ignore
    """Wraps an SOMA DenseNDArray for xarray variable/DataArray support.

    Note: This class does not open or close the SOMA array. The array must be already
    opened for reading to be used with this wrapper.
    """

    def __init__(
        self,
        array: DenseNDArray,
    ):
        self._array = array
        self._dtype: np.typing.DTypeLike = self._array.schema.field(
            "soma_data"
        ).type.to_pandas_dtype()

    def _raw_indexing_method(
        self, key: Tuple[Union[slice, int], ...]
    ) -> np.typing.NDArray[Any]:
        """Returns a numpy array containing data from a requested tuple."""
        assert len(key) == len(self.shape)  # Should be verified by xarray.

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
                raise ValueError("Slice steps are not supported.")
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

    def __getitem__(self, key: ExplicitIndexer) -> np.typing.NDArray[Any]:
        """Returns data from SOMA DenseNDArray using xarray-style indexing."""

        # If any of the slices have steps, convert the key to a vectorized
        # indexer and the slice to a numpy array.
        def has_step(x: Any) -> bool:
            return isinstance(x, slice) and x.step not in (1, None)

        if any(has_step(index) for index in key.tuple):

            key_tuple = tuple(
                (
                    np.arange(index.start, index.stop, index.step)
                    if has_step(index)
                    else index
                )
                for index, dim_size in zip(key.tuple, self.shape)
            )

            key = (
                OuterIndexer(key_tuple)
                if isinstance(key, BasicIndexer) or isinstance(key, OuterIndexer)
                else VectorizedIndexer(key_tuple)
            )

        # Use xarray indexing adapter to further partition indicies.
        return explicit_indexing_adapter(  # type: ignore
            key,
            self.shape,
            IndexingSupport.BASIC,
            self._raw_indexing_method,
        )

    @property
    def dtype(self) -> np.typing.DTypeLike:
        """Numpy dtype of the array data."""
        return self._dtype

    @property
    def shape(self) -> Tuple[int, ...]:
        """Shape of the wrapped SOMA DenseNDArray."""
        return self._array.shape


class DenseNDArrayDatastore(AbstractDataStore):  # type: ignore
    """Xarray datastore for reading a SOMA DenseNDArray.

    This class can be used to read an Xarray ``Dataset`` using the
    ``xarray.Dataset.load_store`` class method. The SOMA DenseNDArray will be opened
    when this class is initialized and closed by the ``close`` method. It will always
    create a ``Dataset`` with exactly one ``DataArray``.

    Args:
        uri: The URI of the :class:`DenseNDarray` to open.
        variable_name: Name to use for the ``DataArray`` holding this
            :class:`DenseNDArray`. Defaults to ``soma_data``.
        dim_names: Name to use for the dimensions. Defaults to the
            :class:`DenseNDArray` dimension names (e.g. ``("soma_dim0", "soma_dim1")``).
        context: The Context value to use when opening the object. Defaults to ``None``.
        platform_config: Platform configuration options specific to this open operation.
            Defaults to ``None``.
    """

    __slots__ = (
        "_array",
        "_attrs",
        "_dim_names",
        "_variable_name",
    )

    def __init__(
        self,
        uri: str,
        *,
        variable_name: str = "soma_data",
        dim_names: Optional[str] = None,
        attrs: Optional[Mapping[str, Any]] = None,
        context: Optional[SOMATileDBContext] = None,
        platform_config: Optional[options.PlatformConfig] = None,
    ):
        """Initialize and open the data store."""
        self._array = DenseNDArray.open(
            uri,
            mode="r",
            context=context,
            platform_config=platform_config,
        )
        self._variable_name = variable_name
        self._dim_names = (
            tuple(f"soma_dim_{enum}" for enum in range(self._array.ndim))
            if dim_names is None
            else dim_names
        )
        self._attrs = Frozen(dict()) if attrs is None else Frozen(attrs)

    def close(self) -> None:
        """Close the data store."""
        self._array.close()

    def load(self) -> Tuple[Frozen[str, Any], Frozen[str, Any]]:
        """Returns a dictionary of ``xarray.Variable`` objects and a dictionary of
        ``xarray.Attribute`` objects.
        """
        variables = Frozen(
            {
                self._variable_name: Variable(
                    self._dim_names,
                    LazilyIndexedArray(DenseNDArrayWrapper(self._array)),
                    attrs={
                        key: val
                        for key, val in self._array.metadata.items()
                        if not key.startswith("soma_")
                    },
                )
            }
        )
        return variables, self._attrs
