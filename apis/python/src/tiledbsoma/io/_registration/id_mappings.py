# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

from typing import Dict, List, Tuple

import anndata as ad
import attrs
import pandas as pd
from typing_extensions import Self

import tiledbsoma
import tiledbsoma.logging


@attrs.define(kw_only=True)
class AxisIDMapping:
    """
    For a single to-be-appended AnnData/H5AD input in SOMA multi-file append-mode ingestion, this
    class tracks the mapping of input-data ``obs`` or ``var`` 0-up offsets to SOMA join ID values
    for the destination SOMA experiment.

    See module-level comments for more information.
    """

    # Tuple not List so this can't be modified by accident when passed into some function somewhere
    data: Tuple[int, ...]

    def is_identity(self) -> bool:
        for i, data in enumerate(self.data):
            if data != i:
                return False
        return True

    def get_shape(self) -> int:
        if len(self.data) == 0:
            return 0
        else:
            return 1 + max(self.data)

    @classmethod
    def identity(cls, n: int) -> Self:
        """This maps 0-up input-file offsets to 0-up soma_joinid values. This is
        important for uns arrays which we never grow on ingest --- rather, we
        sub-nest the entire recursive ``uns`` data structure.
        """
        return cls(data=tuple(range(n)))


@attrs.define(kw_only=True)
class ExperimentIDMapping:
    """
    For a single to-be-appended AnnData/H5AD input in SOMA multi-file append-mode ingestion, this
    class contains an ``ExperimentIDMapping`` for ``obs``, and one ``ExperimentIDMapping`` for
    ``var`` in each measurement.

    See module-level comments for more information.
    """

    obs_axis: AxisIDMapping
    var_axes: Dict[str, AxisIDMapping]

    @classmethod
    def from_isolated_anndata(
        cls,
        adata: ad.AnnData,
        measurement_name: str,
    ) -> Self:
        """Factory method to compute offset-to-SOMA-join-ID mappings for a single input file in
        isolation. This is used when a user is ingesting a single AnnData/H5AD to a single SOMA
        experiment, not in append mode, allowing us to still have the bulk of the ingestor code to
        be non-duplicated between non-append mode and append mode.
        """
        tiledbsoma.logging.logger.info(
            "Registration: registering isolated AnnData object."
        )

        obs_mapping = AxisIDMapping(data=tuple(range(len(adata.obs))))
        var_axes = {}
        var_axes[measurement_name] = AxisIDMapping(data=tuple(range(len(adata.var))))
        if adata.raw is not None:
            var_axes["raw"] = AxisIDMapping(data=tuple(range(len(adata.raw.var))))

        return cls(obs_axis=obs_mapping, var_axes=var_axes)


def get_dataframe_values(df: pd.DataFrame, field_name: str) -> List[str]:
    """Extracts the label values (e.g. cell barcode, gene symbol) from an AnnData/H5AD
    ``obs`` or ``var`` dataframe."""
    if field_name in df:
        values = [str(e) for e in df[field_name]]
    elif df.index.name in (field_name, "index", None):
        values = list(df.index)
    else:
        raise ValueError(f"could not find field name {field_name} in dataframe")

    # Check the values are unique.
    if len(values) != len(set(values)):
        raise ValueError(
            f"non-unique registration values have been provided in field {field_name}"
        )
    return values
