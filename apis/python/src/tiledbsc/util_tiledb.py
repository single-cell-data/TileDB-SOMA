#!/usr/bin/env python

from typing import Optional

import pandas as pd
import tiledb

import tiledbsc


# ----------------------------------------------------------------
def show_soma_schemas(soma_uri: str, ctx: Optional[tiledb.Ctx] = None) -> None:
    """
    Show some summary information about an ingested TileDB Single-Cell Group.  This tool goes a bit beyond `print(tiledb.group.Group(soma_uri))` by also revealing array schema. Additionally, by employing encoded domain-specific knowleldge, it traverses items in the familiar order `X`, `obs`, `var`, etc. rather than using the general-purpose tiledb-group-display function.
    """

    # Tab-completion at the shell can insert a trailing slash; leave it off
    # so we don't show undesired '...//...' in component URIs.
    soma_uri = soma_uri.rstrip("/")
    soma = tiledbsc.SOMA(soma_uri)

    for key in soma.X.keys():
        _show_array_schema(soma.X[key].uri, ctx)
    _show_array_schema(soma.obs.uri, ctx)
    _show_array_schema(soma.var.uri, ctx)

    _show_array_schemas_for_group(soma.obsm.uri, ctx)
    _show_array_schemas_for_group(soma.varm.uri, ctx)
    _show_array_schemas_for_group(soma.obsp.uri, ctx)
    _show_array_schemas_for_group(soma.varp.uri, ctx)

    # Not all groups have raw X data
    if soma.raw.exists():
        for key in soma.raw.X.keys():
            _show_array_schema(soma.raw.X[key].uri, ctx)
            _show_array_schema(soma.raw.var.uri, ctx)


# ----------------------------------------------------------------
def _show_array_schema(uri: str, ctx: Optional[tiledb.Ctx] = None) -> None:
    print("----------------------------------------------------------------")
    print("Array:", uri)
    with tiledb.open(uri, ctx=ctx) as A:
        print(A.schema)


# ----------------------------------------------------------------
def _show_array_schemas_for_group(
    group_uri: str, ctx: Optional[tiledb.Ctx] = None
) -> None:
    with tiledb.Group(group_uri, mode="r", ctx=ctx) as G:
        for element in G:
            if element.type == tiledb.libtiledb.Array:
                _show_array_schema(element.uri, ctx)


# ----------------------------------------------------------------
def show_tiledb_group_array_schemas(uri: str, ctx: Optional[tiledb.Ctx] = None) -> None:
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


# ----------------------------------------------------------------
def list_fragments(array_uri: str) -> None:
    print(f"Listing fragments for array: '{array_uri}'")
    vfs = tiledb.VFS()

    fragments = []
    fi = tiledb.fragment.FragmentInfoList(array_uri=array_uri)

    for f in fi:
        f_dict = {
            "array_schema_name": f.array_schema_name,
            "num": f.num,
            "cell_num": f.cell_num,
            "size": vfs.dir_size(f.uri),
        }

        # parse nonempty domains into separate columns
        for d in range(len(f.nonempty_domain)):
            f_dict[f"d{d}"] = f.nonempty_domain[d]

        fragments.append(f_dict)

    frags_df = pd.DataFrame(fragments)
    print(frags_df)
