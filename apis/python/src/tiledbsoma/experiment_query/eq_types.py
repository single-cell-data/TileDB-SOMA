from typing import Dict, Optional, Sequence

import pyarrow as pa
from typing_extensions import TypedDict


class ExperimentAxisQueryReadArrowResult(TypedDict, total=False):
    """Return type for the ExperimentAxisQuery.read() method."""

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
    total=False,
)
