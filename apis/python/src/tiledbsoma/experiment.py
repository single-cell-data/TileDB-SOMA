from typing import Any, Dict, Optional, Tuple, cast

import somacore
from typing_extensions import Final

from .collection import CollectionBase
from .dataframe import DataFrame
from .measurement import Measurement
from .options import SOMATileDBContext
from .tiledb_object import TileDBObject

_EMPTY_QUERY = somacore.AxisQuery()


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
        # Top-level objects should specify this:
        context: Optional[SOMATileDBContext] = None,
    ):
        """
        Also see the ``TileDBObject`` constructor.
        """
        super().__init__(uri=uri, parent=parent, context=context)

    # Inherited from somacore
    soma_type: Final = "SOMAExperiment"

    def create(self) -> "Experiment":
        """
        Creates the data structure on disk/S3/cloud.
        """
        self._create(self.soma_type)
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
        obs_query: somacore.AxisQuery = _EMPTY_QUERY,
        var_query: somacore.AxisQuery = _EMPTY_QUERY,
    ) -> somacore.ExperimentAxisQuery:
        """
        Create a query on this Experiment. See ``ExperimentAxisQuery`` for more
        information on parameters and usage.
        """
        if not self.exists():
            raise ValueError(f"Experiment {self.uri} does not exist.")
        return somacore.ExperimentAxisQuery(
            # While not technically a somacore.Experiment yet, we implement
            # all the parts that `ExperimentAxisQuery` needs.
            self,  # type: ignore[arg-type]
            measurement_name,
            obs_query=obs_query,
            var_query=var_query,
        )
