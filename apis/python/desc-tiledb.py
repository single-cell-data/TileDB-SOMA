#!/usr/bin/env python

# ================================================================
# This tool goes a bit beyond
#   print(grp._dump(True))
# by also revealing array schema.
# ================================================================

import sys, os
from util_tiledb import show_tiledb_group_array_schemas
from util_tiledb import show_single_cell_group

def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} (one or more TileDB group-directory nmes)", file=sys.stderr)
        sys.exit(1)

    for uri in sys.argv[1:]:
        show_single_cell_group(uri)

if __name__ == "__main__":
    main()
