# See ../tools/peek-soco.py

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

if len(sys.argv) == 1:
    soco_path = "soma-collection"
elif len(sys.argv) == 2:
    soco_path = sys.argv[1]
else:
    print(f"{sys.argv[0]}: need just one soma-collection path.", file=sys.stderr)
    sys.exit(1)

soco = tiledbsoma.SOMACollection(soco_path)

# Interact at the Python prompt now
