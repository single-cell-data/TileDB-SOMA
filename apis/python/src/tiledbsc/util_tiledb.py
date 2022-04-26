#!/usr/bin/env python

import sys, os
import tiledb

def show_single_cell_group(uri):
    """
    Show some summary information about an ingested TileDB Single-Cell Group.
    This goes a bit beyond `print(grp._dump(True))` by also revealing array schema,
    """
    print('================================================================')
    print('X/data:')
    with tiledb.open(uri+'/X/data') as A:
        df = A[:]
        print("keys", [k for k in df.keys()])
        print(df)
        print(A.schema)

    print('----------------------------------------------------------------')
    print('obs:')
    with tiledb.open(uri+'/obs') as A:
        df = A[:]
        print("keys", [k for k in df.keys()])
        print(A.schema)

    print('----------------------------------------------------------------')
    print('var:')
    with tiledb.open(uri+'/var') as A:
        df = A[:]
        print("keys", [k for k in df.keys()])
        print(A.schema)

    for name in ['obsm', 'varm', 'obsp', 'varp']:
        try:
            grp = tiledb.Group(uri+'/'+name, 'r')
            print()
            print('----------------------------------------------------------------')
            print(name, ':', sep='')
            for element in grp:
                with tiledb.open(element.uri) as A:
                    #print(A[:])
                    print("  ", element.uri, A.shape)
                    print(A.schema)
            grp.close()
        except:
            # Not all groups have all four of obsm, obsp, varm, and varp.
            pass

def show_tiledb_group_array_schemas(uri):
    """
    Recursively show array schemas within a TileDB Group. This function is not specific to
    single-cell matrix-API data.
    """
    grp = tiledb.Group(uri, 'r')
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
            with tiledb.open(element.uri) as A:
                print(A.schema)
        else:
            print("Skipping element type", element.type)
    grp.close()
