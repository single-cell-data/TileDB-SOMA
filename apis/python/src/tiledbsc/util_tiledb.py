#!/usr/bin/env python

import sys, os
import tiledb
from typing import Optional

# ================================================================
def show_single_cell_group(soma_uri: str, ctx: Optional[tiledb.Ctx] = None):
    """
    Show some summary information about an ingested TileDB Single-Cell Group.  This tool goes a bit beyond `print(tiledb.group.Group(soma_uri))` by also revealing array schema. Additionally, by employing encoded domain-specific knowleldge, it traverses items in the familiar order `X`, `obs`, `var`, etc. rather than using the general-purpose tiledb-group-display function.
    """

    # Tab-completion at the shell can insert a trailing slash; leave it off
    # so we don't show undesired '...//...' in component URIs.
    soma_uri = soma_uri.rstrip("/")

    __show_array_schema(os.path.join(soma_uri, "X", "data"), ctx)
    __show_array_schema(os.path.join(soma_uri, "obs"), ctx)
    __show_array_schema(os.path.join(soma_uri, "var"), ctx)

    for name in ["obsm", "varm", "obsp", "varp"]:
        __show_array_schemas_for_group(os.path.join(soma_uri, name), ctx)

    # Not all groups have raw X data
    raw_group = None
    raw_group_uri = os.path.join(soma_uri, "raw")
    try:
        raw_group = tiledb.Group(raw_group_uri, mode="r", ctx=ctx)
    except:
        return

    __show_array_schema(os.path.join(raw_group_uri, "X", "data"), ctx)
    __show_array_schema(os.path.join(raw_group_uri, "var"), ctx)
    __show_array_schemas_for_group(os.path.join(raw_group_uri, "varm"), ctx)


# ----------------------------------------------------------------
def __show_array_schema(uri: str, ctx: Optional[tiledb.Ctx] = None):
    print("----------------------------------------------------------------")
    print("Array:", uri)
    with tiledb.open(uri, ctx=ctx) as A:
        print(A.schema)


# ----------------------------------------------------------------
def __show_array_schemas_for_group(group_uri: str, ctx: Optional[tiledb.Ctx] = None):
    with tiledb.Group(group_uri, mode="r", ctx=ctx) as G:
        for element in G:
            if element.type == tiledb.libtiledb.Array:
                __show_array_schema(element.uri, ctx)


# ================================================================
def show_tiledb_group_array_schemas(uri: str, ctx: Optional[tiledb.Ctx] = None):
    """
    Recursively show array schemas within a TileDB Group. This function is not specific to
    single-cell matrix-API data, and won't necessarily traverse items in a familiar
    application-specific order.
    """
    with tiledb.Group(uri, mode="r", ctx=ctx) as G:
        print()
        print("================================================================")
        print(uri)

        for element in G:
            # Note: use `element.type` rather than `isinstance(element, tiledb.group.Group)`
            # since type(element) is `tiledb.object.Object` in all cases.
            if element.type == tiledb.group.Group:
                show_tiledb_group_array_schemas(element.uri)
            elif element.type == tiledb.libtiledb.Array:
                print()
                print(
                    "----------------------------------------------------------------"
                )
                print(element.uri)
                with tiledb.open(element.uri, ctx=ctx) as A:
                    print(A.schema)
            else:
                print("Skipping element type", element.type)
