# Copyright (c) 2024 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2024 TileDB, Inc.
#
# Licensed under the MIT License.

import os

# Temporary for https://github.com/single-cell-data/TileDB-SOMA/issues/2407
_new_shape_feature_flag = os.getenv("SOMA_PY_NEW_SHAPE") is not None


def _new_shape_feature_flag_enabled() -> bool:
    """
    This is temporary only and will be removed once
    https://github.com/single-cell-data/TileDB-SOMA/issues/2407
    is complete.
    """
    return _new_shape_feature_flag
