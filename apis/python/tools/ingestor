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
    parser.add_argument(
        "-q", "--quiet", help="decrease output verbosity", action="store_true"
    )
    parser.add_argument(
        "-n",
        help="All arguments after flags are treated as input paths",
        action="store_true",
    )
    parser.add_argument(
        "-o",
        help="Specify output directory to contain the somas",
        type=str,
        default="./tiledb-data",
    )
    # TODO: add an update option once we figure out tiledb.from_pandas mode=overwrite flag.
    parser.add_argument(
        "--ifexists",
        help="What to do if soma storage already exists. Default: continue.",
        choices=["replace", "abort", "continue"],
        type=str,
        nargs=1,
    )
    parser.add_argument(
        "paths",
        type=str,
        help="One for specified input with default output path, or two to specify input and output paths, or multiple input paths if -n is specified",
        nargs="*",
    )
    args = parser.parse_args()

    if args.ifexists is None:
        args.ifexists = ["continue"]

    outdir = args.o.rstrip("/")

    verbose = not args.quiet

    if args.n:
        if len(args.paths) < 1:
            parser.print_help(file=sys.stderr)
            sys.exit(1)
        for input_path in args.paths:
            # Example 'anndata/pbmc3k_processed.h5ad' -> 'tiledb-data/pbmc3k_processed'
            output_path = os.path.join(
                outdir, os.path.splitext(os.path.basename(input_path))[0]
            )
            ingest_one(input_path, output_path, args.ifexists[0], verbose)
    else:
        if len(args.paths) == 0:
            input_path = "anndata/pbmc-small.h5ad"
            output_path = os.path.join(outdir, "pbmc-small")
            ingest_one(input_path, output_path, args.ifexists[0], verbose)
        elif len(args.paths) == 1:
            input_path = args.paths[0]
            # Example 'anndata/pbmc3k_processed.h5ad' -> 'tiledb-data/pbmc3k_processed'
            output_path = os.path.join(
                outdir, os.path.splitext(os.path.basename(input_path))[0]
            )
            ingest_one(input_path, output_path, args.ifexists[0], verbose)
        elif len(args.paths) == 2:
            input_path = args.paths[0]
            output_path = args.paths[1]
            ingest_one(input_path, output_path, args.ifexists[0], verbose)
        else:
            parser.print_help(file=sys.stderr)
            sys.exit(1)


def ingest_one(input_path: str, output_path: str, ifexists: str, verbose: bool):
    # Check that the input exists.
    if not os.path.exists(input_path):
        # Print this neatly and exit neatly, to avoid a multi-line stack trace otherwise.
        print(f"Input path not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    # Prepare to write the output.
    # This is for local-disk use only -- for S3-backed tiledb://... URIs we should
    # use tiledb.vfs to remove any priors, and/or make use of a tiledb `overwrite` flag.
    parent = os.path.dirname(output_path.rstrip("/"))
    if not os.path.exists(parent):
        os.mkdir(parent)

    if os.path.exists(output_path):
        if ifexists == "continue":
            print(f"Already exists; continuing: {output_path}")
            return
        elif ifexists == "abort":
            print(f"Already exists; aborting: {output_path}")
            sys.exit(1)
        elif ifexists == "replace":
            print(f"Already exists; replacing: {output_path}")
            shutil.rmtree(output_path)  # Overwrite
        else:
            raise Exception(
                "Internal coding error in --ifexists handling.", ifexists, "<"
            )

    soma = tiledbsc.SOMA(uri=output_path, verbose=verbose)

    # Do the ingest into TileDB.
    soma.from_h5ad(input_path)
    if not verbose:
        print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
