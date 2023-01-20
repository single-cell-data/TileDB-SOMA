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


def count_obs(exp: tiledbsoma.Experiment, attr_name: str) -> None:
    print(
        exp.obs.read(column_names=[attr_name])
        .concat()
        .to_pandas()
        .groupby(attr_name)
        .size()
        .sort_values()
    )


if len(sys.argv) == 1:
    input_path = "tiledbsoma-data/pbmc-small"
    # input_path = 'tiledbsoma-data/pbmc3k_processed'
elif len(sys.argv) == 2:
    input_path = sys.argv[1]
else:
    print(f"{sys.argv[0]}: need just one Experiment path.", file=sys.stderr)
    sys.exit(1)

exp = tiledbsoma.Experiment(input_path)
if not exp.exists():
    print("Does not exist yet:", input_path)

# Interact at the Python prompt now
