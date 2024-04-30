# Copyright (c) 2024 TileDB, Inc,
#
# Licensed under the MIT License.
"""Implementation of a SOMA image collections."""

from typing import Any, Optional, Sequence

import numpy as np
import pyarrow as pa
import somacore
from somacore import coordinates, images, options

from . import _tdb_handles
from ._collection import CollectionBase
from ._dense_nd_array import DenseNDArray
from ._soma_object import AnySOMAObject


class Image2D(  # type: ignore[misc]  # __eq__ false positive
    CollectionBase[AnySOMAObject],
    images.Image2D[DenseNDArray, AnySOMAObject],
):
    """TODO: Add documentaiton for Image2D

    Lifecycle:
        Experimental.
    """

    __slots__ = ()
    _wrapper_type = _tdb_handles.Image2DWrapper

    def add_new_level(
        self, key: str, *, scale: Sequence[np.float64], **kwargs: Any
    ) -> DenseNDArray:
        raise NotImplementedError()

    @property
    def level_count(self) -> int:
        raise NotImplementedError()

    def level_properties(self, level: int) -> images.Image2D.LevelProperties:
        raise NotImplementedError()

    def read_level(
        self,
        level: int,
        coords: options.DenseNDCoords = (),
        *,
        transform: Optional[coordinates.CoordinateTransform] = None,
        result_order: options.ResultOrderStr = somacore.ResultOrder.ROW_MAJOR,
        platform_config: Optional[options.PlatformConfig] = None,
    ) -> pa.Tensor:
        """TODO: Add read_image_level documentation"""
        raise NotImplementedError()
