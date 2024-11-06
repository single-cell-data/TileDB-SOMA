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

DENSE_ARRAYS_CAN_HAVE_CURRENT_DOMAIN = clib.embedded_version_triple() >= (2, 27, 0)

# Temporary for # https://github.com/single-cell-data/TileDB-SOMA/issues/2407:
# this allows testing dense + current domain on the same machine without
# having to switch core builds (or switch machines).
if os.getenv("SOMA_IGNORE_CORE_2_27") is not None:
    DENSE_ARRAYS_CAN_HAVE_CURRENT_DOMAIN = False
