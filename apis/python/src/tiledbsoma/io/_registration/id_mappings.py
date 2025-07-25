# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.
from __future__ import annotations

from typing import cast

import anndata as ad
import attrs
import numpy as np
import numpy.typing as npt
import pandas as pd
from typing_extensions import Self


@attrs.define(kw_only=True, frozen=True)
class AxisIDMapping:
    """For a single to-be-appended AnnData/H5AD input in SOMA multi-file append-mode ingestion, this
    class tracks the mapping of input-data ``obs`` or ``var`` 0-up offsets to SOMA join ID values
    for the destination SOMA experiment.

    Private class
    """

    data: npt.NDArray[np.int64]

    def __attrs_post_init__(self) -> None:
        self.data.setflags(write=False)

    def get_shape(self) -> int:
        if len(self.data) == 0:
            return 0
        else:
            return int(self.data.max() + 1)

    def is_identity(self) -> bool:
        # fast rejection first
        if self.get_shape() != len(self.data) or self.data[0] != 0:
            return False

        return np.array_equal(self.data, np.arange(0, len(self.data)))

    @classmethod
    def identity(cls, n: int) -> Self:
        """This maps 0-up input-file offsets to 0-up soma_joinid values. This is
        important for uns arrays which we never grow on ingest --- rather, we
        sub-nest the entire recursive ``uns`` data structure.
        """
        return cls(data=np.arange(n, dtype=np.int64))


@attrs.define(kw_only=True, frozen=True)
class ExperimentIDMapping:
    """For a single to-be-appended AnnData/H5AD input in SOMA multi-file append-mode ingestion, this
    class contains an ``AxisIDMapping`` for ``obs``, and one ``AxisIDMapping`` for
    ``var`` in each measurement.

    Private class
    """

    obs_axis: AxisIDMapping
    var_axes: dict[str, AxisIDMapping]

    @classmethod
    def from_anndata(cls, adata: ad.AnnData, *, measurement_name: str = "RNA") -> Self:
        """Create a new ID mapping from an AnnData.

        This is useful for creating a new Experiment from a single AnnData.
        """
        obs_axis = AxisIDMapping.identity(len(adata.obs))
        var_axes = {measurement_name: AxisIDMapping.identity(len(adata.var))}
        if adata.raw is not None:
            var_axes["raw"] = AxisIDMapping.identity(len(adata.raw.var))
        return cls(obs_axis=obs_axis, var_axes=var_axes)


def get_dataframe_values(df: pd.DataFrame, field_name: str) -> pd.Series:
    """Extracts the label values (e.g. cell barcode, gene symbol) from an AnnData/H5AD
    ``obs`` or ``var`` dataframe.
    """
    if field_name in df:
        values = cast(pd.Series, df[field_name].astype(str))
    elif df.index.name in (field_name, "index", None):
        values = cast(pd.Series, df.index.to_series().astype(str))
    else:
        raise ValueError(f"Could not find field name {field_name} in dataframe.")

    # Check the values are unique.
    if not values.is_unique:
        raise ValueError(
            f"Non-unique registration values have been provided in field {field_name}."
        )
    return values
