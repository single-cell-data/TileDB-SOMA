#!/bin/bash

# See ../repl/peek-exp

import sys

import anndata
import numpy
import pandas
import scipy  # noqa: F401
import tiledb  # noqa: F401

import tiledbsoma
import tiledbsoma.io

# module aliases
ad = anndata
np = numpy
pd = pandas


def count_obs(exp: tiledbsoma.SOMAExperiment, attr_name: str) -> None:
    print(
        exp.obs.to_dataframe(attrs=[attr_name]).groupby(attr_name).size().sort_values()
    )


if len(sys.argv) == 1:
    input_path = "tiledbsoma-data/pbmc-small"
    # input_path = 'tiledbsoma-data/pbmc3k_processed'
elif len(sys.argv) == 2:
    input_path = sys.argv[1]
else:
    print(f"{sys.argv[0]}: need just one SOMAExperiment path.", file=sys.stderr)
    sys.exit(1)

cfg = tiledb.Config()
cfg["py.init_buffer_bytes"] = 4 * 1024**3
ctx = tiledb.Ctx(cfg)

exp = tiledbsoma.SOMAExperiment(input_path, ctx=ctx)
if not exp.exists():
    print("Does not exist yet:", input_path)

# Interact at the Python prompt now
