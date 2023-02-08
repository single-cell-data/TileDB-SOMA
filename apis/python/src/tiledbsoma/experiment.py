from somacore import experiment

from .collection import Collection, CollectionBase
from .dataframe import DataFrame
from .measurement import Measurement
from .tiledb_object import AnyTileDBObject


class Experiment(
    CollectionBase[AnyTileDBObject],
    experiment.Experiment[  # type: ignore[type-var]
        DataFrame,
        Collection[Measurement],
        AnyTileDBObject,
    ],
):
    """
    ``obs``: Primary annotations on the observation axis. The contents of the
             ``soma_joinid`` column define the observation index domain,
             AKA ``obs_id``. All observations for the Experiment must be
             defined in this dataframe.

    ``ms``: A collection of named measurements.

    [lifecycle: experimental]
    """

    __slots__ = ()

    _subclass_constrained_soma_types = {
        "obs": ("SOMADataFrame",),
        "ms": ("SOMACollection",),
    }
