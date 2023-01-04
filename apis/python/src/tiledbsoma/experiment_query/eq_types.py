"""
Types global to this module
"""
from typing import Dict, Optional, Sequence

import pyarrow as pa
from typing_extensions import TypedDict

# Sadly, you can't define a generic TypedDict....


class ExperimentQueryReadArrowResult(TypedDict, total=False):
    obs: pa.Table
    var: pa.Table
    X: pa.Table
    X_layers: Dict[str, pa.Table]


AxisColumnNames = TypedDict(
    "AxisColumnNames",
    {
        "obs": Optional[Sequence[str]],  # None is all
        "var": Optional[Sequence[str]],
    },
)
