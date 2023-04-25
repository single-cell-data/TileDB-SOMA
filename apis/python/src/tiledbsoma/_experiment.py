# Copyright (c) 2021-2023 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2023 TileDB, Inc.
#
# Licensed under the MIT License.

"""Implementation of a SOMA Experiment.
"""

from typing import Any

from somacore import experiment

from ._collection import Collection, CollectionBase
from ._dataframe import DataFrame
from ._measurement import Measurement
from ._tdb_handles import Wrapper
from ._tiledb_object import AnyTileDBObject


class Experiment(
    CollectionBase[AnyTileDBObject],
    experiment.Experiment[  # type: ignore[type-var]
        DataFrame,
        Collection[Measurement],
        AnyTileDBObject,
    ],
):
    """An :class:`Experiment` represents a single-cell experiment. It is
    a container class that has observations and measurements.

    Attributes:
        obs (DataFrame):
            Primary annotations on the observation axis. The contents of the
            ``soma_joinid`` column define the observation index domain,
            AKA ``obs_id``. All observations for the Experiment must be
            defined in this dataframe.
        ms (Collection):
            A collection of named measurements.

    Lifecycle:
        Experimental.
    """

    __slots__ = ()

    _subclass_constrained_soma_types = {
        "obs": ("SOMADataFrame",),
        "ms": ("SOMACollection",),
    }

    @classmethod
    def _set_create_metadata(cls, handle: Wrapper[Any]) -> None:
        # Root SOMA objects include a `dataset_type` entry to allow the
        # TileDB Cloud UI to detect that they are SOMA datasets.
        handle.metadata["dataset_type"] = "soma"
        return super()._set_create_metadata(handle)
