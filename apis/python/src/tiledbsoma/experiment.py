from typing import Any, Dict, Optional, Tuple, cast

import tiledb
from typing_extensions import Final

from .collection import CollectionBase
from .dataframe import DataFrame
from .experiment_query import AxisQuery, ExperimentAxisQuery
from .measurement import Measurement
from .tiledb_object import TileDBObject
from .tiledb_platform_config import TileDBPlatformConfig


class Experiment(CollectionBase[TileDBObject]):
    """
    ``obs``: Primary annotations on the observation axis. The contents of the
             ``soma_joinid`` column define the observation index domain,
             AKA ``obs_id``. All observations for the Experiment must be
             defined in this dataframe.

    ``ms``: A collection of named measurements.
    """

    _subclass_constrained_soma_types: Dict[str, Tuple[str, ...]] = {
        "obs": ("SOMADataFrame",),
        "ms": ("SOMACollection",),
    }

    def __init__(
        self,
        uri: str,
        *,
        # Non-top-level objects can have a parent to propagate context, depth, etc.
        parent: Optional[CollectionBase[Any]] = None,
        # Top-level objects should specify these:
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

    soma_type: Final = "SOMAExperiment"

    def create(self) -> "Experiment":
        """
        Creates the data structure on disk/S3/cloud.
        """
        super().create()
        return self

    @property
    def obs(self) -> DataFrame:
        """
        Primary annotations on the observation axis. The contents of the
        ``soma_joinid`` column define the observation index domain,
        AKA ``obs_id``. All observations for the Experiment must be
        defined in this dataframe.
        """
        return cast(DataFrame, self["obs"])

    @property
    def ms(self) -> CollectionBase[Measurement]:
        """
        A collection of named measurements.
        """
        return cast(CollectionBase[Measurement], self["ms"])

    def axis_query(
        self,
        measurement_name: str,
        *,
        obs_query: Optional[AxisQuery] = None,
        var_query: Optional[AxisQuery] = None,
    ) -> ExperimentAxisQuery:
        """
        Create a query on this Experiment. See ``ExperimentAxisQuery`` for more
        information on parameters and usage.
        """
        return ExperimentAxisQuery(
            self, measurement_name, obs_query=obs_query, var_query=var_query
        )
