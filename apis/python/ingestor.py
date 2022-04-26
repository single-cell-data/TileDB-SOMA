#!/usr/bin/env python

# ================================================================
# A simple driver for ingestion of anndata to a TileDB group.
#
# * Invoke this with no arguments to process:
#   o canonical input:  ./anndata/pbmc3k_processed.h5ad
#   o canonical output: ./tiledb-data/pbmc3k_processed/
#
# * Invoke this with one argument /path/to/some/somename.h5ad:
#   o Output will be ./tiledb-data/somename
#
# * Invoke this with two arguments to specify input anndata HDF5 file
#   and output TileDB group.
#
# Nominal immediate-term support is to local disk, although output to tiledb:/...
# URIs will be supported.
#
# Note this removes and recreates the destination TileDB group on each invocation.
# ================================================================

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

scdataset = tiledbsc.SCGroup(output_path, verbose=True)
scdataset.from_h5ad(input_path)
