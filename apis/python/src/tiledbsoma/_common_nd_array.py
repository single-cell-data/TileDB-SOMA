# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""Common code shared by both NDArray implementations."""

from __future__ import annotations

from collections.abc import Sequence
from typing import cast

import pyarrow as pa
from typing_extensions import Self

from ._core_options import PlatformConfig
from ._soma_array import SOMAArray
from ._soma_context import SOMAContext
from ._types import OpenTimestamp, StatusAndReason
from .options._soma_tiledb_context import SOMATileDBContext
from .options._tiledb_create_write_options import TileDBCreateOptions


class NDArray(SOMAArray):
    """Abstract base for the common behaviors of both kinds of NDArray."""

    __slots__ = ()

    @classmethod
    def create(
        cls,
        uri: str,
        *,
        type: pa.DataType,
        shape: Sequence[int | None],
        platform_config: PlatformConfig | None = None,
        context: SOMAContext | SOMATileDBContext | None = None,
        tiledb_timestamp: OpenTimestamp | None = None,
    ) -> Self:
        """Creates a SOMA ``NDArray`` at the given URI.

        Args:
            type:
                The Arrow type to be stored in the NDArray.
                If the type is unsupported, an error will be raised.
            shape:
                The maximum capacity of each dimension, including room
                for any intended future appends, as a sequence.  E.g. ``(100, 10)``.
                All lengths must be in the positive int64 range.

                For :class:`SparseNDArray` only, if a slot is None, then the minimum
                possible range will be used.  This makes a :class:`SparseNDArray`
                growable using ``resize``.
            platform_config:
                Platform-specific options used to create this array.
                This may be provided as settings in a dictionary, with options
                located in the ``{'tiledb': {'create': ...}}`` key,
                or as a :class:`~tiledbsoma.TileDBCreateOptions` object.
            context: If provided, the :class:`SOMAContext` to use when creating and opening this collection. If not,
                provide the default context will be used and possibly initialized.
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
            TileDBError:
                If unable to create the underlying object.

        Lifecycle:
            Maturing.
        """
        raise NotImplementedError("must be implemented by child class.")

    def resize(self, newshape: Sequence[int | None], check_only: bool = False) -> StatusAndReason:
        """Increases the shape of the array as specfied. Raises an error if the new
        shape is less than the current shape in any dimension. Raises an error if
        the new shape exceeds maxshape in any dimension. Raises an error if the
        array doesn't already have a shape: in that case please call
        tiledbsoma_upgrade_shape. If ``check_only`` is ``True``, returns
        whether the operation would succeed if attempted, and a reason why it
        would not.

        Lifecycle:
            Maturing.
        """
        if check_only:
            return cast("StatusAndReason", self._handle.can_resize(newshape))
        self._handle.resize(newshape)
        return (True, "")

    def tiledbsoma_upgrade_shape(self, newshape: Sequence[int | None], check_only: bool = False) -> StatusAndReason:
        """Allows the array to have a resizeable shape as described in the TileDB-SOMA
        1.15 release notes.  Raises an error if the new shape exceeds maxshape in
        any dimension. Raises an error if the array already has a shape.
        """
        if check_only:
            return cast("StatusAndReason", self._handle.tiledbsoma_can_upgrade_shape(newshape))
        self._handle.tiledbsoma_upgrade_shape(newshape)
        return (True, "")

    @property
    def shape(self) -> tuple[int, ...]:
        """Returns capacity of each dimension, always a list of length ``ndim``.
        This will not necessarily match the bounds of occupied cells within the array.
        Rather, it is the bounds outside of which no data may be read or written.

        Lifecycle:
            Maturing.
        """
        return cast("tuple[int, ...]", tuple(self._handle.shape))

    @property
    def maxshape(self) -> tuple[int, ...]:
        """Returns the maximum resizable capacity of each dimension, always a list of length
        ``ndim``.  This will not necessarily match the bounds of occupied cells within the array.
        It is the upper limit for ``resize`` on the array.

        Lifecycle:
            Maturing.
        """
        return cast("tuple[int, ...]", tuple(self._handle.maxshape))

    @property
    def ndim(self) -> int:
        """The number of dimensions in this array.

        Lifecycle:
            Maturing.
        """
        return len(self.shape)

    @property
    def tiledbsoma_has_upgraded_shape(self) -> bool:
        """Returns true if the array has the upgraded resizeable shape feature
        from TileDB-SOMA 1.15: the array was created with this support, or it has
        had ``.tiledbsoma_upgrade_shape`` applied to it.

        Lifecycle:
            Maturing.
        """
        return cast("bool", self._handle.tiledbsoma_has_upgraded_shape)

    @property
    def type(self) -> pa.DataType:
        """Returns the data type of the array's soma_data attribute.

        This is a read-only property that provides direct access to the Arrow data type
        of the array's data, equivalent to ``self.schema.field("soma_data").type``.

        Lifecycle:
            Maturing.

        Returns:
            The Arrow data type of the array's data.
        """
        return self.schema.field("soma_data").type

    @classmethod
    def _dim_capacity_and_extent(
        cls,
        dim_name: str,
        dim_shape: int | None,
        ndim: int,
        create_options: TileDBCreateOptions,
    ) -> tuple[int, int]:
        raise NotImplementedError("must be implemented by child class.")
