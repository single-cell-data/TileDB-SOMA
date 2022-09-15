# See ../tools/peek-soma.py

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


def count_obs(soma: tiledbsoma.SOMA, attr_name: str) -> None:
    print(soma.obs.df(attrs=[attr_name]).groupby(attr_name).size().sort_values())


def count_var(soma: tiledbsoma.SOMA, attr_name: str) -> None:
    print(soma.var.df(attrs=[attr_name]).groupby(attr_name).size().sort_values())


if len(sys.argv) == 1:
    input_path = "tiledb-data/pbmc-small"
    # input_path = 'tiledb-data/pbmc3k_processed'
elif len(sys.argv) == 2:
    input_path = sys.argv[1]
else:
    print(f"{sys.argv[0]}: need just one soma path.", file=sys.stderr)
    sys.exit(1)

soma = tiledbsoma.SOMA(input_path)
if not soma.exists():
    print("Does not exist yet:", input_path)

# Interact at the Python prompt now
