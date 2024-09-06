# Copyright (c) 2024 TileDB, Inc.
#
# Licensed under the MIT License.

import json
from typing import Dict, List, Literal, Union, cast

from .pytiledbsoma import tiledbsoma_stats_string

ParsedStats = List[Dict[Literal["counters", "timers"], Dict[str, Union[float, int]]]]


def tiledbsoma_stats_json() -> str:
    """Returns tiledbsoma stats as a JSON string"""
    # cast is needed for pybind11 things
    return cast(str, tiledbsoma_stats_string())


def tiledbsoma_stats_as_py() -> ParsedStats:
    """Returns tiledbsoma stats as a Python dict"""
    # cast is needed for pybind11 things
    return cast(ParsedStats, json.loads(tiledbsoma_stats_string()))
