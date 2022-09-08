# See ../repl/peek-soma

import sys

import anndata
import numpy
import pandas
import scipy  # noqa: F401
import tiledb  # noqa: F401

import tiledbsc
import tiledbsc.io

# module aliases
ad = anndata
np = numpy
pd = pandas


def count_obs(soma: tiledbsc.SOMA, attr_name: str) -> None:
    print(soma.obs.df(attrs=[attr_name]).groupby(attr_name).size().sort_values())


def count_var(soma: tiledbsc.SOMA, attr_name: str) -> None:
    print(soma.var.df(attrs=[attr_name]).groupby(attr_name).size().sort_values())


if len(sys.argv) == 1:
    input_path = "tiledb-data/pbmc-small"
    # input_path = 'tiledb-data/pbmc3k_processed'
elif len(sys.argv) == 2:
    input_path = sys.argv[1]
else:
    print(f"{sys.argv[0]}: need just one soma path.", file=sys.stderr)
    sys.exit(1)

cfg = tiledb.Config()
cfg["py.init_buffer_bytes"] = 4 * 1024**3
ctx = tiledb.Ctx(cfg)

soma = tiledbsc.SOMA(input_path, ctx=ctx)
if not soma.exists():
    print("Does not exist yet:", input_path)

# Interact at the Python prompt now
