from typing import Dict, Optional, Tuple

import somacore
import somacore.experiment as scexp

from .collection import CollectionBase
from .options import SOMATileDBContext
from .tiledb_object import TileDBObject

_EMPTY_QUERY = somacore.AxisQuery()


class Experiment(scexp.Experiment[TileDBObject], CollectionBase[TileDBObject]):
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
        context: Optional[SOMATileDBContext] = None,
    ):
        """
        Also see the ``TileDBObject`` constructor.
        """
        super().__init__(uri=uri, context=context)

    # Inherited from somacore
    # soma_type: Final = "SOMAExperiment"

    def _legacy_create(self) -> "Experiment":
        """
        Creates the data structure on disk/S3/cloud.
        """
        self._create(self.soma_type)
        return self

    def axis_query(
        self,
        measurement_name: str,
        *,
        obs_query: somacore.AxisQuery = _EMPTY_QUERY,
        var_query: somacore.AxisQuery = _EMPTY_QUERY,
    ) -> somacore.ExperimentAxisQuery["Experiment"]:  # type: ignore[type-var]
        """
        Create a query on this Experiment. See ``ExperimentAxisQuery`` for more
        information on parameters and usage.
        """
        if not self.exists():
            raise ValueError(f"Experiment {self.uri} does not exist.")
        return super().axis_query(
            measurement_name, obs_query=obs_query, var_query=var_query
        )
