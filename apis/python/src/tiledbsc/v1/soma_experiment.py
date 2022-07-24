from typing import Any, Dict, Optional

import tiledb

from .soma_collection import SOMACollection
from .soma_dataframe import SOMADataFrame
from .tiledb_object import TileDBObject
from .tiledb_platform_config import TileDBPlatformConfig


class SOMAExperiment(SOMACollection):
    """
    `obs`: Primary annotations on the _observation_ axis. The contents of the `soma_rowid` pseudo-column define
    the _observation_ index domain, aka `obsid`. All observations for the SOMAExperiment _must_ be
    defined in this dataframe.

    `ms`: A collection of named measurements.
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
        Also see the :class:`TileDBObject` constructor.
        """
        super().__init__(
            uri=uri,
            name=name,
            parent=parent,
            tiledb_platform_config=tiledb_platform_config,
            ctx=ctx,
        )
        self._constructors = {
            "obs": SOMADataFrame,
            "ms": SOMACollection,
        }
        self._cached_members = {}

    def create(self) -> None:
        """
        Creates the data structure on disk/S3/cloud.
        """
        super().create()

    def __getattr__(self, name: str) -> Any:  # TODO: union type
        """
        Implements `experiment.obs` and `experiment.ms`.
        """
        if name in self._constructors:
            if name not in self._cached_members:
                child_uri = self._get_child_uri(name)
                self._cached_members[name] = self._constructors[name](
                    uri=child_uri, name=name, parent=self
                )
            return self._cached_members[name]

        else:
            # Unlike __getattribute__ this is _only_ called when the member isn't otherwise
            # resolvable. So raising here is the right thing to do.
            raise AttributeError(f"unrecognized attribute: {name}")
