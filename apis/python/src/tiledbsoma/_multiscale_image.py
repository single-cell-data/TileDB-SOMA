#
# Copyright (c) 2024 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2024 TileDB, Inc.
#
# Licensed under the MIT License.
"""
Implementation of a SOMA MultiscaleImage.
"""

import warnings
from typing import Any, Optional, Sequence, Tuple, Union

import pyarrow as pa
import somacore
from somacore import CoordinateSpace, CoordinateTransform, ScaleTransform, options
from typing_extensions import Self

from ._constants import SPATIAL_DISCLAIMER
from ._dense_nd_array import DenseNDArray
from ._soma_object import AnySOMAObject
from .options import SOMATileDBContext


class ImageProperties:
    """Properties for a single resolution level in a multiscale image.

    Lifecycle:
        Experimental.
    """

    @property
    def name(self) -> str:
        """The key for the image.

        Lifecycle:
            Experimental
        """
        raise NotImplementedError()

    @property
    def shape(self) -> Tuple[int, ...]:
        """Size of each axis of the image.

        Lifecycle:
            Experimental.
        """
        raise NotImplementedError()


class MultiscaleImage(somacore.MultiscaleImage[DenseNDArray, AnySOMAObject]):
    """A multiscale image with an extendable number of resolution levels.
    """A multiscale image represented as a collection of images at multiple resolution levels.

    Each level of the multiscale image must have the following consistent properties:

    * **Number of Channels**: All levels must have the same number of channels.
    * **Axis Order**: The order of axes (e.g., channels, height, width) must be consistent across levels.

    Lifecycle:
        Experimental.
    """

    __slots__ = ()

    # Lifecycle

    @classmethod
    def create(
        cls,
        uri: str,
        *,
        type: pa.DataType,
        reference_level_shape: Sequence[int],
        axis_names: Sequence[str] = ("c", "y", "x"),
        axis_types: Sequence[str] = ("channel", "height", "width"),
        platform_config: Optional[options.PlatformConfig] = None,
        context: Optional[SOMATileDBContext] = None,
    ) -> Self:
        """Creates a new collection of this type at the given URI.

        Args:
            uri: The URI where the collection will be created.
            reference_level_shape: The shape of the reference level for the multiscale
                image. In most cases, this corresponds to the size of the image
                at ``level=0``.
            axis_names: The names of the axes of the image.
            axis_types: The types of the axes of the image. Must be the same length as
                ``axis_names``. Valid types are: ``channel``, ``height``, ``width``,
                and ``depth``.

        Returns:
            The newly created collection, opened for writing.

        Lifecycle:
            Experimental.
        """
        warnings.warn(SPATIAL_DISCLAIMER)
        raise NotImplementedError()

    def add_new_level(
        self,
        key: str,
        *,
        uri: Optional[str] = None,
        shape: Sequence[int],
        **kwargs: Any,
    ) -> DenseNDArray:
        """Add a new level in the multi-scale image.

        Parameters are as in :meth:`DenseNDArray.create`. The provided shape will
        be used to compute the scale between images and must correspond to the image
        size for the entire image.

        Lifecycle:
            Experimental.
        """
        raise NotImplementedError()

    # Data operations

    def read_spatial_region(
        self,
        level: Union[int, str],
        region: Optional[options.SpatialRegion] = None,
        *,
        channel_coords: options.DenseCoord = None,
        region_transform: Optional[CoordinateTransform] = None,
        region_coord_space: Optional[CoordinateSpace] = None,
        result_order: options.ResultOrderStr = options.ResultOrder.ROW_MAJOR,
        platform_config: Optional[options.PlatformConfig] = None,
    ) -> somacore.SpatialRead[pa.Tensor]:
        """Reads a user-defined region of space into a :class:`SpatialRead` with data
        in either an Arrow tensor or table.

        Reads the bounding box of the input region from the requested image level. This
        will return a :class:`SpatialRead` with the image data stored as a
        :class:`pa.Tensor`.

        Args:
            level: The image level to read the data from. May use index of the level
                or the image name.
            region: The region to query. May be a box in the form
                [x_min, y_min, x_max, y_max] (for 2D images), a box in the form
                [x_min, y_min, z_min, x_max, y_max, z_max] (for 3D images), or
                a shapely Geometry.
            channel_coords: An optional slice that defines the channel coordinates
                to read.
            region_transform: An optional coordinate transform that provides the
                transformation from the provided region to the reference level of this
                image. Defaults to ``None``.
            region_coord_space: An optional coordinate space for the region being read.
                The axis names must match the input axis names of the transform.
                Defaults to ``None``, coordinate space will be inferred from transform.
            result_order: the order to return results, specified as a
                :class:`~options.ResultOrder` or its string value.

        Returns:
            The data bounding the requested region as a :class:`SpatialRead` with
            :class:`pa.Tensor` data.

        Lifecycle:
            Experimental.
        """
        raise NotImplementedError()

    # Metadata operations

    @property
    def axis_names(self) -> Tuple[str, ...]:
        """The name of the image axes.

        Lifecycle:
            Experimental.
        """
        raise NotImplementedError()

    @property
    def coordinate_space(self) -> Optional[CoordinateSpace]:
        """Coordinate space for this multiscale image.

        Lifecycle:
            Experimental.
        """
        raise NotImplementedError()

    @coordinate_space.setter
    def coordinate_space(self, value: CoordinateSpace) -> None:
        """Coordinate space for this multiscale image.

        Lifecycle:
            Experimental.
        """
        raise NotImplementedError()

    def get_transform_from_level(self, level: Union[int, str]) -> ScaleTransform:
        """Returns the transformation from user requested level to image reference
        level.

        Lifecycle:
            Experimental.
        """
        raise NotImplementedError()

    def get_transform_to_level(self, level: Union[int, str]) -> ScaleTransform:
        """Returns the transformation from the image reference level to the user
        requested level.

        Lifecycle:
            Experimental.
        """
        raise NotImplementedError()

    @property
    def image_type(self) -> str:
        """The order of the axes as stored in the data model.

        Lifecycle:
            Experimental.
        """
        raise NotImplementedError()

    @property
    def level_count(self) -> int:
        """The number of image levels stored in the ``MultiscaleImage``.

        Lifecycle:
            Experimental.
        """
        raise NotImplementedError()

    def level_properties(self, level: Union[int, str]) -> somacore.ImageProperties:
        """The properties of an image at the specified level.

        Lifecycle:
            Experimental.
        """
        raise NotImplementedError()

    @property
    def reference_level(self) -> Optional[int]:
        """The index of image level that is used as a reference level.

        This will return ``None`` if no current image level matches the size of the
        reference level.

        Lifecycle:
            Experimental.
        """
        raise NotImplementedError()

    @property
    def reference_level_properties(self) -> "ImageProperties":
        """The image properties of the reference level.

        Lifecycle:
            Experimental.
        """
        raise NotImplementedError()
