from typing import Dict, Tuple

from somacore import experiment

from . import collection
from .tiledb_object import AnyTileDBObject


class Experiment(  # type: ignore[misc]
    experiment.Experiment[AnyTileDBObject],
    collection.AnyTileDBCollection,
):
    """
    ``obs``: Primary annotations on the observation axis. The contents of the
             ``soma_joinid`` column define the observation index domain,
             AKA ``obs_id``. All observations for the Experiment must be
             defined in this dataframe.

    ``ms``: A collection of named measurements.

    [lifecycle: experimental]
    """

    _subclass_constrained_soma_types: Dict[str, Tuple[str, ...]] = {
        "obs": ("SOMADataFrame",),
        "ms": ("SOMACollection",),
    }
