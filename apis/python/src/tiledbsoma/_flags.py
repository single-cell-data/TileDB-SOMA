# Copyright (c) 2024 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2024 TileDB, Inc.
#
# Licensed under the MIT License.

import os

import tiledbsoma.pytiledbsoma as clib

# This is temporary for
# https://github.com/single-cell-data/TileDB-SOMA/issues/2407.  It will be
# removed once https://github.com/single-cell-data/TileDB-SOMA/issues/2407 is
# complete.

NEW_SHAPE_FEATURE_FLAG_ENABLED = os.getenv("SOMA_PY_NEW_SHAPE") != "false"

DENSE_ARRAYS_CAN_HAVE_CURRENT_DOMAIN = clib.embedded_version_triple() >= (2, 27, 0)
