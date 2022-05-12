#!/usr/bin/env python

# ================================================================
# This tool goes a bit beyond
#   print(tiledb.group.Group('tiledb-data/pbmc3k_processed'))
# by also revealing array schema.
# ================================================================

import sys, os
import tiledbsc

def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} (one or more TileDB group-directory nmes)", file=sys.stderr)
        sys.exit(1)

    for uri in sys.argv[1:]:
        tiledbsc.util_tiledb.show_single_cell_group(uri)

if __name__ == "__main__":
    main()
