#!/usr/bin/env python

import pathlib

TEST_DIR = pathlib.Path(__file__).parent
SOMA_URI = f"{TEST_DIR}/../../../data/soco/pbmc3k_processed"

if not pathlib.Path(SOMA_URI).exists():
    raise RuntimeError("Please run `make data` in the repo base directory")
