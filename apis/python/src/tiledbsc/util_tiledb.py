#!/usr/bin/env python

import sys, os
import tiledb
from typing import Optional

def show_single_cell_group(uri: str, ctx: Optional[tiledb.Ctx] = None):
    """
    Show some summary information about an ingested TileDB Single-Cell Group.
    This tool goes a bit beyond
      `print(tiledb.group.Group('tiledb-data/pbmc3k_processed')._dump(True))`
    by also revealing array schema.
    """

    print('================================================================')
    print('X/data:')
    with tiledb.open(os.path.join(uri, 'X', 'data'), ctx=ctx) as A:
        df = A[:]
        print("keys", list(df.keys()))
        print(df)
        print(A.schema)

    # Not all groups have raw X data
    try:
        with tiledb.open(uri+'/X/raw', ctx=ctx) as A:
            print('X/raw:')
            df = A[:]
            print("keys", list(df.keys()))
            print(df)
            print(A.schema)
    except:
        pass

    print('----------------------------------------------------------------')
    print('obs:')
    with tiledb.open(os.path.join(uri, 'obs'), ctx=ctx) as A:
        df = A[:]
        print("keys", list(df.keys()))
        print(A.schema)

    print('----------------------------------------------------------------')
    print('var:')
    with tiledb.open(os.path.join(uri, 'var'), ctx=ctx) as A:
        df = A[:]
        print("keys", list(df.keys()))
        print(A.schema)

    for name in ['obsm', 'varm', 'obsp', 'varp']:
        # Not all groups have all four of obsm, obsp, varm, and varp.
        grp = None
        try:
            grp = tiledb.Group(os.path.join(uri, name), mode='r', ctx=ctx)
        except:
            pass

        if grp != None:
            print()
            print('----------------------------------------------------------------')
            print(name, ':', sep='')
            for element in grp:
                with tiledb.open(element.uri, ctx=ctx) as A:
                    print(element.uri)
                    print(A.schema)
            grp.close()
            pass

def show_tiledb_group_array_schemas(uri: str, ctx: Optional[tiledb.Ctx] = None):
    """
    Recursively show array schemas within a TileDB Group. This function is not specific to
    single-cell matrix-API data.
    """
    grp = tiledb.Group(uri, mode='r', ctx=ctx)
    print()
    print('================================================================')
    print(uri)

    for element in grp:
        # Note: use `element.type` rather than `isinstance(element, tiledb.group.Group)`
        # since type(element) is `tiledb.object.Object` in all cases.
        if element.type == tiledb.group.Group:
            show_tiledb_group_array_schemas(element.uri)
        elif element.type == tiledb.libtiledb.Array:
            print()
            print('----------------------------------------------------------------')
            print(element.uri)
            with tiledb.open(element.uri, ctx=ctx) as A:
                print(A.schema)
        else:
            print("Skipping element type", element.type)
    grp.close()
