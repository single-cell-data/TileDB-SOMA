from typing import Any, Dict, Literal, Optional, Tuple, cast

import tiledb

from .collection import CollectionBase
from .dataframe import DataFrame
from .measurement import Measurement
from .tiledb_object import TileDBObject
from .tiledb_platform_config import TileDBPlatformConfig


class Experiment(CollectionBase[TileDBObject]):
    """
    ``obs``: Primary annotations on the observation axis. The contents of the
             ``soma_joinid`` column define the observation index domain,
             aka ``obs_id``. All observations for the Experiment must be
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

    @property
    def soma_type(self) -> Literal["SOMAExperiment"]:
        return "SOMAExperiment"

    def create(self) -> "Experiment":
        """
        Creates the data structure on disk/S3/cloud.
        """
        super().create()
        return self

    @property
    def obs(self) -> DataFrame:
        return cast(DataFrame, self["obs"])

    @property
    def ms(self) -> CollectionBase[Measurement]:
        return cast(CollectionBase[Measurement], self["ms"])
