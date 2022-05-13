#!/usr/bin/env python

# ================================================================
# A simple driver for ingestion of anndata to a TileDB group.
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
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Ingest soma data from anndata/h5ad into TileDB group structure"
    )
    parser.add_argument("-q", "--quiet", help="decrease output verbosity", action="store_true")
    parser.add_argument(
        "paths",
        type=str,
        help="One for specified input with default output path, or two to specify input and output paths, or multiple input paths if -n is specified",
        nargs='*'
    )
    parser.add_argument("-n", help="All arguments after flags are treated as input paths", action="store_true")
    args = parser.parse_args()

    verbose = not args.quiet

    if args.n:
        if len(args.paths) < 1:
            parser.print_help(file=sys.stderr)
            sys.exit(1)
        for input_path in args.paths:
            # Example 'anndata/pbmc3k_processed.h5ad' -> 'tiledb-data/pbmc3k_processed'
            output_path = 'tiledb-data/' + os.path.splitext(os.path.basename(input_path))[0]
            ingest_one(input_path, output_path, verbose)
    else:
        if len(args.paths) == 0:
            input_path  = 'anndata/pbmc-small.h5ad'
            output_path = 'tiledb-data/pbmc-small'
            ingest_one(input_path, output_path, verbose)
        elif len(args.paths) == 1:
            input_path  = args.paths[0]
            # Example 'anndata/pbmc3k_processed.h5ad' -> 'tiledb-data/pbmc3k_processed'
            output_path = 'tiledb-data/' + os.path.splitext(os.path.basename(input_path))[0]
            ingest_one(input_path, output_path, verbose)
        elif len(args.paths) == 2:
            input_path  = args.paths[0]
            output_path = args.paths[1]
            ingest_one(input_path, output_path, verbose)
        else:
            parser.print_help(file=sys.stderr)
            sys.exit(1)


def ingest_one(input_path: str, output_path: str, verbose: bool):
    # Check that the input exists.
    if not os.path.exists(input_path):
        # Print this neatly and exit neatly, to avoid a multi-line stack trace otherwise.
        print(f"Input path not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    # Prepare to write the output.
    # This is for local-disk use only -- for S3-backed tiledb://... URIs we should
    # use tiledb.vfs to remove any priors, and/or make use of a tiledb `overwrite` flag.
    parent = os.path.dirname(output_path.rstrip('/'))
    if not os.path.exists(parent):
        os.mkdir(parent)
    # TODO: make this an option, with default off. Currently we need an overwrite-flag for
    # tiledb.from-numpy in UnsArray in order to allow users to do update-in-place.
    if os.path.exists(output_path):
        shutil.rmtree(output_path) # Overwrite

    soma = tiledbsc.SOMA(uri=output_path, verbose=verbose)

    # Do the ingest into TileDB.
    soma.from_h5ad(input_path)
    if not verbose:
        print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
