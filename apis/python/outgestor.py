#!/usr/bin/env python

# ================================================================
# A simple driver for outgestion of a TileDB soma group to anndata.
#
# * Invoke this with one argument /path/to/some/somagroup:
#   o Output will be ./anndata-readback/somagroup.h5ad
#
# * Invoke this with two arguments to specify input TileDB group and
#   output anndata HDF5 file.
#
# Nominal immediate-term support is to local disk, although input from tiledb:/...
# URIs will be supported.
# ================================================================

import tiledb
import tiledbsc

import anndata as ad

import pandas as pd
import scipy
import numpy as np

import sys, os, shutil
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Outgest soma data from TileDB group structure to anndata/h5ad"
   )
    parser.add_argument("-q", "--quiet", help="decrease output verbosity", action="store_true")
    parser.add_argument(
        "paths",
        type=str,
        help="One for specified input with default output path, or two to specify input and output paths",
        nargs='*'
    )
    args = parser.parse_args()

    if len(args.paths) == 0:
        input_path  = 'tiledb-data/pbmc-small'
        output_path = 'anndata-readback/pbmc-small.h5ad'
    elif len(args.paths) == 1:
        # Strip trailing slashes so basename will behave correctly
        input_path  = args.paths[0].rstrip('/')
        # Example 'tiledb-data/pbmc3k_processed' -> 'anndata-readcbak/pbmc3k_processed.h5ad'
        output_path = 'anndata-readback/' + os.path.basename(input_path) + '.h5ad'
    elif len(args.paths) == 2:
        input_path  = args.paths[0]
        output_path = args.paths[1]
    else:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    if not os.path.exists('anndata-readback'):
        os.mkdir('anndata-readback')

    if not os.path.exists(input_path):
        # Print this neatly and exit neatly, to avoid a multi-line stack trace otherwise.
        print(f"Input path not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    # This is for local-disk use only -- for S3-backed tiledb://... URIs we should
    # use tiledb.vfs to remove any priors, and/or make use of a tiledb `overwrite` flag.
    if not os.path.exists('anndata-readback'):
        os.mkdir('anndata-readback')

    verbose = not args.quiet
    soma = tiledbsc.SOMA(input_path, verbose=verbose)
    soma.to_h5ad(output_path)

    if not verbose:
        print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
