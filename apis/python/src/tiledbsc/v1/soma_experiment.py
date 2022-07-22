from typing import Optional

from .soma_collection import SOMACollection


class SOMAExperiment(SOMACollection):
    """
    TBD
    """

    # TODO

    # `obs`
    # `SOMADataFrame`
    # Primary annotations on the _observation_ axis. The contents of the `__rowid` pseudo-column define
    # the _observation_ index domain, aka `obsid`. All observations for the SOMAExperiment _must_ be
    # defined in this dataframe.

    # `ms`
    # `SOMACollection[string, SOMAMeasurement]`
    # A collection of named measurements.

    def __init__(
        self,
        uri: str,
        *,
        name: Optional[str] = None,
        parent: Optional[SOMACollection] = None,
    ):
        """
        See also the :class:`TileDBOject` constructor.
        """
        super().__init__(uri=uri, name=name, parent=parent)

        # TODO: try to make these lazily instantiated, rather than here in the constructor.  This is
        # important for cloud performance, as name-to-URI resolution always involves one or more
        # HTTP requests.

        # See comments in _get_child_uris
        # child_uris = self._get_child_uris(["obs", "ms"])

        # obs_uri = child_uris["obs"]
        # ms_uri = child_uris["ms"]

        # self.obs = SOMADatatFrame(uri=obs_uri, name="obs", parent=self)
        # self.ms = SOMACollection(uri=ms_uri, name="ms", parent=self)
