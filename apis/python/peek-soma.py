#!/usr/bin/env python -i

# Invoke this with, for example,
#
#   python -i peek-soma.py tiledb-data/pbmc3k_processed
#
# -- then you can inspect the soma object.

import tiledb
import tiledbsc
import sys, os

import anndata
import anndata as ad # so we can type it either way

import pandas
import pandas as pd # so we can type it either way
import numpy
import numpy  as np # so we can type it either way
import scipy

if len(sys.argv) == 1:
    input_path = 'tiledb-data/pbmc-small'
    #input_path = 'tiledb-data/pbmc3k_processed'
elif len(sys.argv) == 2:
    input_path  = sys.argv[1]
else:
    print(f"{sys.argv[0]}: need just one soma path.", file=sys.stderr)
    sys.exit(1)

soma = tiledbsc.SOMA(input_path)

# Interact at the Python prompt now
