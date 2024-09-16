# Copyright (c) 2021-2023 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2023 TileDB, Inc.
#
# Licensed under the MIT License.

"""Common code shared by both NDArray implementations."""

from typing import Optional, Sequence, Tuple, Union, cast

import pyarrow as pa
import somacore
from somacore import options
from typing_extensions import Self

from ._soma_array import SOMAArray
from ._types import OpenTimestamp
from .options._soma_tiledb_context import (
    SOMATileDBContext,
)
from .options._tiledb_create_write_options import TileDBCreateOptions


class NDArray(SOMAArray, somacore.NDArray):
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
                All lengths must be in the positive int64 range, or ``None``.  It's
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
            tiledbsoma.AlreadyExistsError:
                If the underlying object already exists at the given URI.
            tiledbsoma.NotCreateableError:
                If the URI is malformed for a particular storage backend.
            TileDBError:
                If unable to create the underlying object.

        Lifecycle:
            Maturing.
        """
        raise NotImplementedError("must be implemented by child class.")

    @property
    def shape(self) -> Tuple[int, ...]:
        """Returns capacity of each dimension, always a list of length ``ndim``.
        This will not necessarily match the bounds of occupied cells within the array.
        Rather, it is the bounds outside of which no data may be written.

        Lifecycle:
            Maturing.
        """
        return cast(Tuple[int, ...], tuple(self._handle.shape))

    @property
    def maxshape(self) -> Tuple[int, ...]:
        """Returns the maximum resizable capacity of each dimension, always a list of length
        ``ndim``.  This will not necessarily match the bounds of occupied cells within the array.
        It is the upper limit for ``resize`` on the array.

        Lifecycle:
            Maturing.
        """
        return cast(Tuple[int, ...], tuple(self._handle.maxshape))

    @property
    def tiledbsoma_has_upgraded_shape(self) -> bool:
        """Returns true if the array has the upgraded resizeable shape feature
        from TileDB-SOMA 1.15: the array was created with this support, or it has
        had ``.tiledbsoma_upgrade_shape`` applied to it.

        Lifecycle:
            Maturing.
        """
        return self._handle.tiledbsoma_has_upgraded_shape

    @classmethod
    def _dim_capacity_and_extent(
        cls,
        dim_name: str,
        dim_shape: Optional[int],
        create_options: TileDBCreateOptions,
    ) -> Tuple[int, int]:
        raise NotImplementedError("must be implemented by child class.")
