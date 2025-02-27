#!/usr/bin/env python

import os

TEST_DIR = os.path.dirname(__file__)
SOMA_URI = f"{TEST_DIR}/../../../data/soco/pbmc3k_processed"

if not os.path.exists(SOMA_URI):
    raise RuntimeError("Please run `make data` in the repo base directory")
