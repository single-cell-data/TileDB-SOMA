#!/usr/bin/env python

# ================================================================
# This is an anndata-describer that goes a bit beyond what h5ls does for us.
#
# See also:
#
# * `brew install hdf5`
# * `h5ls -r anndata/pbmc3k_processed.h5ad`
# * `h5ls -vr anndata/pbmc3k_processed.h5ad`
#
# Please see comments in util_ann.py.
# ================================================================

import sys
from util_ann import describe_ann_file

def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} (one or more ANN data-file names)", file=sys.stderr)
        sys.exit(1)

    for input_path in sys.argv[1:]:
        describe_ann_file(input_path)

if __name__ == "__main__":
    main()
