from typing import Any, Optional

import tiledb

from .soma_collection import SOMACollection
from .soma_dataframe import SOMADataFrame
from .tiledb_platform_config import TileDBPlatformConfig


class SOMAExperiment(SOMACollection):
    """
    TBD
    """

    """
    Primary annotations on the _observation_ axis. The contents of the `soma_rowid` pseudo-column define
    the _observation_ index domain, aka `obsid`. All observations for the SOMAExperiment _must_ be
    defined in this dataframe.
    """
    _cached_obs: Optional[SOMADataFrame]

    """
    A collection of named measurements.
    """
    _cached_ms: Optional[SOMACollection]  # of SOMAMeasurement

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
        Also see the :class:`TileDBObject` constructor.
        """
        super().__init__(
            uri=uri,
            name=name,
            parent=parent,
            tiledb_platform_config=tiledb_platform_config,
            ctx=ctx,
        )
        self._cached_obs = None
        self._cached_ms = None

    def create(self) -> None:
        """
        Creates the data structure on disk/S3/cloud.
        """
        super().create()

    def __getattr__(self, name: str) -> Any:  # TODO: union type
        """
        TODO: COMMENT
        """
        if name == "obs":
            if self._cached_obs is None:
                child_uri = self._get_child_uri("obs")
                self._cached_obs = SOMADataFrame(uri=child_uri, name="obs", parent=self)
            return self._cached_obs
        elif name == "ms":
            if self._cached_ms is None:
                child_uri = self._get_child_uri("ms")
                self._cached_ms = SOMACollection(uri=child_uri, name="ms", parent=self)
            return self._cached_ms
        else:
            # Unlike __getattribute__ this is _only_ called when the member isn't otherwise
            # resolvable. So raising here is the right thing to do.
            raise AttributeError(f"unrecognized attribute: {name}")
