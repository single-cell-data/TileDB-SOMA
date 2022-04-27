#!/usr/bin/env python

# Invoke this with, for example,
#
#   python -i interactor.py ~/scdata/anndata/pbmc-small.h5ad tempdir
#
# -- then you can inspect the anndata object.

import tiledbsc
import sys, os, shutil

input_path  = 'anndata/pbmc3k_processed.h5ad'
output_path = 'tiledb-data/pbmc3k_processed'
if len(sys.argv) == 3:
    input_path  = sys.argv[1]
    output_path = sys.argv[2]
if len(sys.argv) == 2:
    input_path  = sys.argv[1]
    # Example 'anndata/pbmc3k_processed.h5ad' -> 'tiledb-data/pbmc3k_processed'
    output_path = 'tiledb-data/' + os.path.splitext(os.path.basename(input_path))[0]

# This is for local-disk use only -- for S3-backed tiledb://... URIs we should
# use tiledb.vfs to remove any priors, and/or make use of a tiledb `overwrite` flag.
if not os.path.exists('tiledb-data'):
    os.mkdir('tiledb-data')
if os.path.exists(output_path):
    shutil.rmtree(output_path) # Overwrite

scgroup = tiledbsc.SCGroup(output_path, verbose=True)

anndata = scgroup.read_h5ad(input_path)
anndata = scgroup.decategoricalize(anndata)
# Interact at the prompt now:
# * anndata.X
# * anndata.obs.keys()
# * anndata.obs.head()
# * anndata.obs.index
# * anndata.obs.values
# * etc
