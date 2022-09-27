from typing import Any, Dict, Optional

import tiledb

from .soma_collection import SOMACollection
from .soma_dataframe import SOMADataFrame
from .soma_measurement import SOMAMeasurement
from .tiledb_object import TileDBObject
from .tiledb_platform_config import TileDBPlatformConfig


class SOMAExperiment(SOMACollection):
    """
    ``obs``: Primary annotations on the observation axis. The contents of the ``soma_rowid`` pseudo-column define the observation index domain, aka ``obsid``. All observations for the SOMAExperiment must be defined in this dataframe.

    ``ms``: A collection of named measurements.
    """

    _constructors: Dict[str, Any]
    _cached_members: Dict[str, TileDBObject]

    def __init__(
        self,
        uri: str,
        *,
        name: Optional[str] = None,
        # Non-top-level objects can have a parent to propagate context, depth, etc.
        parent: Optional[SOMACollection] = None,
        # Top-level objects should specify these:
        tiledb_platform_config: Optional[TileDBPlatformConfig] = None,
        ctx: Optional[tiledb.Ctx] = None,
    ):
        """
        Also see the ``TileDBObject`` constructor.
        """
        super().__init__(
            uri=uri,
            name=name,
            parent=parent,
            tiledb_platform_config=tiledb_platform_config,
            ctx=ctx,
        )
        self._constructors = {
            # TODO: union-type of SOMADataFrame and SOMAIndexedDataFrame
            "obs": SOMADataFrame,
            "ms": SOMACollection,
        }
        self._cached_members = {}

    def create(self) -> None:
        """
        Creates the data structure on disk/S3/cloud.
        """
        super().create()

    @property
    def obs(self) -> Any:
        return self["obs"]

    @property
    def ms(self) -> Any:
        return self["ms"]

    def __getitem__(self, name: str) -> Any:  # TODO: union type
        """
        Implements ``experiment.obs`` and ``experiment.ms``.
        """
        if name in self._constructors:
            if name not in self._cached_members:
                child_uri = self._get_child_uri(name)
                self._cached_members[name] = self._constructors[name](
                    uri=child_uri, name=name, parent=self
                )
            return self._cached_members[name]

        # otherwise let generic collection handle it.
        super().__getitem__(name)

    def constrain(self) -> None:
        """
        Checks constraints on the ``SOMAExperiment``. Raises an exception if any is violated.
        """
        # TODO: find a good spot to call this from.

        for element in self.ms:
            if not isinstance(element, SOMAMeasurement):
                raise Exception(
                    f"element {element.name} of {self.type}.ms should be SOMAMeasurement; got {element.__class__.__name__}"
                )
            element.constrain()
