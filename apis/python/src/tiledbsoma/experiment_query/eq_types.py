from typing import Dict, Optional, Sequence

import pyarrow as pa
from typing_extensions import TypedDict


class ExperimentAxisQueryReadArrowResult(TypedDict, total=False):
    """Return type for the ExperimentAxisQuery.read() method"""

    obs: pa.Table
    """Experiment.obs query slice, as an Arrow Table"""
    var: pa.Table
    """Experiment.ms[...].var query slice, as an Arrow Table"""
    X: pa.Table
    """Experiment.ms[...].X[...] query slice, as an Arrow Table"""
    X_layers: Dict[str, pa.Table]
    """Any additional X layers requested, as Arrow Table(s)"""


class AxisColumnNames(TypedDict, total=False):
    """Specify column names for the ExperimentAxisQuery API read operations"""

    obs: Optional[Sequence[str]]  # None is all columns
    var: Optional[Sequence[str]]  # None is all columns
