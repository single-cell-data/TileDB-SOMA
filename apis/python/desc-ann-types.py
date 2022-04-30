#!/usr/bin/env python

# ================================================================
# This is an anndata-describer that goes a bit beyond what h5ls does for us.
# Similar to desc-ann.py but only shows essential data-type information.
#
# Please see comments in util_ann.py.
# ================================================================

import sys
import tiledbsc

def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} (one or more ANN data-file names)", file=sys.stderr)
        sys.exit(1)

    for input_path in sys.argv[1:]:
        tiledbsc.util_ann.describe_ann_file(input_path, True)

if __name__ == "__main__":
    main()
