from dataclasses import dataclass
from typing import Dict, List

import anndata as ad
import pandas as pd
from typing_extensions import Self

import tiledbsoma
import tiledbsoma.logging


@dataclass
class AxisIDMapping:
    """TODO: docstring"""

    data: List[int]

    @classmethod
    def identity(cls, n: int) -> Self:
        """TODO: docstring"""
        return cls(list(range(n)))


@dataclass
class ExperimentIDMapping:
    """TODO: docstring"""

    obs_axis: AxisIDMapping
    var_axes: Dict[str, AxisIDMapping]

    @classmethod
    def from_isolated_anndata(
        cls,
        adata: ad.AnnData,
        measurement_name: str,
    ) -> Self:
        """TODO: docstring"""
        tiledbsoma.logging.logger.info(
            "Registration: registering isolated AnnData object."
        )

        obs_mapping = AxisIDMapping(list(range(len(adata.obs))))
        var_axes = {}
        var_axes[measurement_name] = AxisIDMapping(list(range(len(adata.var))))
        if adata.raw is not None:
            var_axes["raw"] = AxisIDMapping(list(range(len(adata.raw.var))))

        return cls(obs_axis=obs_mapping, var_axes=var_axes)


def get_dataframe_values(df: pd.DataFrame, field_name: str) -> List[str]:
    """TODO: docstring"""
    if field_name in df:
        return [str(e) for e in df[field_name]]
    if field_name == df.index.name:
        return list(df.index)
    if df.index.name is None:
        return list(df.index)
    raise ValueError(f"could not find field name {field_name} in dataframe")
