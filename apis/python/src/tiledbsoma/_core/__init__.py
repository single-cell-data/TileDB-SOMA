# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.
"""Module base for the Python reference specification of SOMA.

Types will be defined in their own modules and then imported here for a single
unified namespace.
"""

# TODO: pyarrow >= 14.0.1 doesn't play well with some other PyPI packages
# on Mac OS: https://github.com/apache/arrow/issues/42154
# Remove this once we can pin to recent pyarrow.
import pyarrow_hotfix  # noqa: F401

from .base import SOMAObject
from .collection import Collection
from .coordinates import (
    AffineTransform,
    Axis,
    CoordinateSpace,
    CoordinateTransform,
    IdentityTransform,
    ScaleTransform,
    UniformScaleTransform,
)
from .data import DataFrame, DenseNDArray, NDArray, ReadIter, SparseNDArray, SparseRead
from .experiment import Experiment
from .measurement import Measurement
from .options import BatchSize, IOfN, ResultOrder
from .query import AxisColumnNames, AxisQuery, ExperimentAxisQuery
from .scene import Scene
from .spatial import GeometryDataFrame, MultiscaleImage, PointCloudDataFrame, SpatialRead
from .types import ContextBase

__all__ = (
    "AffineTransform",
    "Axis",
    "AxisColumnNames",
    "AxisQuery",
    "BatchSize",
    "Collection",
    "ContextBase",
    "CoordinateSpace",
    "CoordinateTransform",
    "DataFrame",
    "DenseNDArray",
    "Experiment",
    "ExperimentAxisQuery",
    "GeometryDataFrame",
    "IOfN",
    "IdentityTransform",
    "Measurement",
    "MultiscaleImage",
    "NDArray",
    "PointCloudDataFrame",
    "ReadIter",
    "ResultOrder",
    "SOMAObject",
    "ScaleTransform",
    "Scene",
    "SparseNDArray",
    "SparseRead",
    "SpatialRead",
    "UniformScaleTransform",
)
