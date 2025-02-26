# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""Implementation of a SOMA Experiment.
"""

from __future__ import annotations

import functools

from somacore import experiment, query

from . import _tdb_handles
from ._collection import Collection, CollectionBase
from ._dataframe import DataFrame
from ._indexer import IntIndexer
from ._measurement import Measurement
from ._query import ExperimentAxisQuery
from ._scene import Scene
from ._soma_object import AnySOMAObject


class Experiment(  # type: ignore[misc]  # `SOMAObject.__eq__`, `SOMAGroup.set` false positives
    CollectionBase[AnySOMAObject],
    experiment.Experiment[  # type: ignore[type-var]
        DataFrame,
        Collection[Measurement],
        Collection[Scene],
        AnySOMAObject,
    ],
):
    """A collection subtype that combines observations and measurements
    from an individual experiment.

    In single cell biology, this can represent multiple modes of measurement
    across a single collection of cells (i.e., a "multimodal dataset").
    Within an experiment, a set of measurements on a single set of variables
    (i.e., features) is represented as a :class:`Measurement`.

    Attributes:
        obs (DataFrame):
            Primary annotations on the observation axis. The contents of the
            ``soma_joinid`` column define the observation index domain,
            AKA ``obs_id``. All observations for the Experiment must be
            defined in this dataframe.
        ms (Collection):
            A collection of named measurements.
        spatial (Collection):
            A collection of spatial scenes.

    Example:
        >>> import tiledbsoma
        >>> with tiledbsoma.open("/path/to/experiment") as exp:
        ...     # While users can interact directly with an Experiment's fields:
        ...     obs_df = exp.obs
        ...
        ...     # the primary use case is to run queries on the experiment data.
        ...     q = exp.axis_query(
        ...         "mtdna",
        ...         obs_query=tiledbsoma.AxisQuery(value_filter="tissue == 'lung'"),
        ...         var_query=tiledbsoma.AxisQuery(coords=(slice(50, 100),)),
        ...     )
        ...     query_obs = q.obs().concat().to_pandas()
        ...     query_var = q.var().concat().to_pandas()

    Lifecycle:
        Maturing.
    """

    __slots__ = ()
    _wrapper_type = _tdb_handles.ExperimentWrapper

    _subclass_constrained_soma_types = {
        "obs": ("SOMADataFrame",),
        "ms": ("SOMACollection",),
        "spatial": ("SOMACollection",),
        "obs_spatial_presence": ("SOMADataFrame",),
    }

    def axis_query(  # type: ignore
        self,
        measurement_name: str,
        *,
        obs_query: query.AxisQuery | None = None,
        var_query: query.AxisQuery | None = None,
    ) -> ExperimentAxisQuery:
        """Creates an axis query over this experiment.
        Lifecycle: Maturing.
        """
        return ExperimentAxisQuery(
            self,
            measurement_name,
            obs_query=obs_query or query.AxisQuery(),
            var_query=var_query or query.AxisQuery(),
            index_factory=functools.partial(
                IntIndexer,
                context=self.context,
            ),
        )
