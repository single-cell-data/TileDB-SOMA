"""
Types global to this module
"""
from typing import Dict, Optional, Sequence

import pandas as pd
import pyarrow as pa
from typing_extensions import TypedDict

# Sadly, you can't define a generic TypedDict....


class ExperimentQueryReadArrowResult(TypedDict, total=False):
    obs: pa.Table
    var: pa.Table
    X: pa.Table
    X_layers: Dict[str, pa.Table]


class ExperimentQueryReadPandasResult(TypedDict, total=False):
    obs: pd.DataFrame
    var: pd.DataFrame
    X: pd.DataFrame
    X_layers: Dict[str, pd.DataFrame]


AxisColumnNames = TypedDict(
    "AxisColumnNames",
    {
        "obs": Optional[Sequence[str]],  # None is all
        "var": Optional[Sequence[str]],
    },
)
