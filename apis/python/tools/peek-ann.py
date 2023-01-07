# See ../tools/peek-ann.py

import sys

import anndata
import numpy
import pandas
import scipy  # noqa: F401
import tiledb  # noqa: F401

# module aliases
ad = anndata
np = numpy
pd = pandas

if len(sys.argv) == 1:
    input_path = "anndata/pbmc-small.h5ad"
elif len(sys.argv) == 2:
    input_path = sys.argv[1]
else:
    print(f"{sys.argv[0]}: need just one .h5ad file name.", file=sys.stderr)
    sys.exit(1)

ann = anndata.read_h5ad(input_path, "r")


# Interact at the prompt now:
# * ann.X
# * ann.obs.keys()
# * ann.obs.head()
# * ann.obs.index
# * ann.obs.values
# * etc
