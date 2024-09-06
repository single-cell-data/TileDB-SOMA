from typing import Any, Optional, Sequence, Tuple

import pyarrow as pa
import somacore
from somacore import CoordinateSpace, options
from typing_extensions import Self

from ._constants import SOMA_GEOMETRY, SOMA_JOINID
from ._read_iters import TableReadIter
from ._spatial_dataframe import SpatialDataFrame

_UNBATCHED = options.BatchSize()


class GeometryDataFrame(SpatialDataFrame, somacore.GeometryDataFrame):
    """A multi-column table of geometries with spatial indexing and a user-defined
    schema.

    Lifecycle: experimental
    """

    __slots__ = ()

    @classmethod
    def create(
        cls,
        uri: str,
        *,
        schema: pa.Schema,
        index_column_names: Sequence[str] = (SOMA_JOINID, SOMA_GEOMETRY),
        axis_names: Sequence[str] = ("x", "y"),
        domain: Optional[Sequence[Optional[Tuple[Any, Any]]]] = None,
        platform_config: Optional[options.PlatformConfig] = None,
        context: Optional[Any] = None,
    ) -> Self:
        """TODO: Add docstring"""
        raise NotImplementedError()

    # Data operations

    def count(self) -> int:
        raise NotImplementedError()

    def read(
        self,
        coords: options.SparseDFCoords = (),
        column_names: Optional[Sequence[str]] = None,
        *,
        result_order: options.ResultOrderStr = options.ResultOrder.AUTO,
        value_filter: Optional[str] = None,
        batch_size: options.BatchSize = _UNBATCHED,
        partitions: Optional[options.ReadPartitions] = None,
        platform_config: Optional[options.PlatformConfig] = None,
    ) -> TableReadIter:
        raise NotImplementedError()

    def write(
        self, values: pa.Table, platform_config: Optional[options.PlatformConfig] = None
    ) -> Self:
        raise NotImplementedError()

    # Metadata operations

    @property
    def axis_names(self) -> Tuple[str, ...]:
        raise NotImplementedError()

    @property
    def coordinate_space(self) -> CoordinateSpace:
        raise NotImplementedError()
