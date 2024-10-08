# Copyright (c) 2024 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2024 TileDB, Inc.
#
# Licensed under the MIT License.

"""
Provides support for the new shape feature in TileDB-SOMA 1.15, including the
ability to process all dataframes/arrays contained within a TileDB-SOMA
Experiment. Please also see
https://github.com/single-cell-data/TileDB-SOMA/issues/2407.  """

from typing import Any, Dict, Optional, Tuple, TypedDict

import tiledbsoma


class SizingArgs(TypedDict):
    """Convenience type-alias for kwargs passed to experiment-level
    upgrade/resize functions.
    """

    nobs: Optional[int]
    nvars: Optional[Dict[str, int]]
    ms_name: Optional[str]
    coll_name: Optional[str]
    verbose: bool
    check_only: bool
    context: Optional[tiledbsoma.SOMATileDBContext]


def show_experiment_shapes(
    uri: str,
    *,
    context: Optional[tiledbsoma.SOMATileDBContext] = None,
) -> bool:
    """For each dataframe/array contained within the SOMA ``Experiment`` pointed
    to by the given URI, shows the deprecated ``used_shape`` (for N-D arrays) or
    the ``count`` (for dataframes), along with the ``shape`` and ``maxshape``
    (for arrays) or ``domain`` and ``maxdomain`` (for dataframes).

    XXX TO DO: each parameter and its description.
    """
    args: SizingArgs = dict(
        nobs=None,
        nvars=None,
        ms_name=None,
        coll_name=None,
        verbose=True,
        check_only=True,
        context=context,
    )

    ok = _treewalk(
        uri,
        visitor=_leaf_visitor_show_shapes,
        args=args,
    )
    return ok


def upgrade_experiment_shapes(
    uri: str,
    *,
    verbose: bool = False,
    check_only: bool = False,
    context: Optional[tiledbsoma.SOMATileDBContext] = None,
) -> bool:
    """For each dataframe contained within the SOMA ``Experiment`` pointed to by
    the given URI, sets the ``domain`` to match the dataframe's current
    ``count``.  For each N-D array, sets the ``shape`` to match the array's
    current ``used_shape``. If ``verbose`` is set to ``True``, an activity log
    is printed. If ``check_only`` is true, only does a dry run and reports any
    reasons the upgrade would fail.

    XXX TO DO: each parameter and its description.
    """
    args: SizingArgs = dict(
        nobs=None,
        nvars=None,
        ms_name=None,
        coll_name=None,
        verbose=verbose,
        check_only=check_only,
        context=context,
    )
    ok = _treewalk(
        uri,
        visitor=_leaf_visitor_upgrade,
        args=args,
    )
    return ok


def resize_experiment(
    uri: str,
    *,
    nobs: int,
    nvars: Dict[str, int],
    verbose: bool = False,
    check_only: bool = False,
    context: Optional[tiledbsoma.SOMATileDBContext] = None,
) -> bool:
    """For each dataframe contained within the SOMA ``Experiment`` pointed to by
    the given URI, resizes the ``domain`` for the ``soma_joinid index column
    (if it is an indexed column) to match the desired new value.

    The desired new value may be the same size as at present, or bigger, but not
    exceeding ``maxdomain`` for the ``soma_joinid`` index column.

    For each N-D array contained within the SOMA ``Experiment``, resizes the ``shape``
    to match the desired values.

    * For ``X`` arrays, the resize is to new ``nobs`` x the measurement's new ``nvar``.
    * For ``obsm`` arrays, the resize is to new ``nobs`` x the array's existing ``soma_dim_1`` shape.
    * For ``obsp`` arrays, the resize is to new ``nobs`` x that same ``nobs``.
    * For ``varm`` arrays, the resize is to new ``nvar`` x the array's existing ``soma_dim_1`` shape.
    * For ``varp`` arrays, the resize is to new ``nvar`` x that same ``nvar``.

    In all cases, the desired new ``shape`` value may be the same size on the
    given dimension, or bigger, but not exceeding ``maxshape`` for the given
    dimension.

    If any array has not been upgraded, then its resize will fail.

    XXX TO DO: each parameter and its description.
    """
    args: SizingArgs = dict(
        nobs=nobs,
        nvars=nvars,
        ms_name=None,
        coll_name=None,
        verbose=verbose,
        check_only=check_only,
        context=context,
    )

    # XXX early check on keys of exp.ms vs keys of nvars

    ok = _treewalk(
        uri,
        visitor=_leaf_visitor_resize,
        args=args,
    )
    return ok


def _treewalk(
    uri: str,
    *,
    node_name: Optional[str] = None,
    visitor: Any,
    args: SizingArgs,
) -> bool:
    retval = True
    with tiledbsoma.open(uri) as item:

        if isinstance(item, tiledbsoma.Experiment):
            if "obs" in item:
                ok = _treewalk(
                    item["obs"].uri, node_name="obs", visitor=visitor, args=args
                )
                retval = retval and ok
            if "ms" in item:
                ok = _treewalk(
                    item["ms"].uri, node_name="ms", visitor=visitor, args=args
                )
                retval = retval and ok

        elif isinstance(item, tiledbsoma.Measurement):
            if "var" in item:
                ok = _treewalk(
                    item["var"].uri, node_name="var", visitor=visitor, args=args
                )
                retval = retval and ok

            for coll_name in ["X", "obsm", "obsp", "varm", "varp"]:
                if coll_name in item:

                    args = dict(
                        nobs=args["nobs"],
                        nvars=args["nvars"],
                        ms_name=args["ms_name"],
                        coll_name=coll_name,
                        verbose=args["verbose"],
                        check_only=args["check_only"],
                        context=args["context"],
                    )

                    ok = _treewalk(
                        item[coll_name].uri,
                        node_name=coll_name,
                        visitor=visitor,
                        args=args,
                    )
                    retval = retval and ok

        elif isinstance(item, tiledbsoma.Collection):

            for key in item:

                if node_name == "ms":
                    args = dict(
                        nobs=args["nobs"],
                        nvars=args["nvars"],
                        ms_name=key,
                        coll_name=args["coll_name"],
                        verbose=args["verbose"],
                        check_only=args["check_only"],
                        context=args["context"],
                    )

                ok = _treewalk(item[key].uri, node_name=key, visitor=visitor, args=args)
                retval = retval and ok

        else:
            ok = visitor(item, node_name=node_name, args=args)
            retval = retval and ok

    return retval


def _leaf_visitor_show_shapes(
    item: Any,
    *,
    node_name: str,
    args: SizingArgs,
) -> bool:
    retval = True
    if isinstance(item, tiledbsoma.DataFrame):
        _print_leaf_node_banner("DataFrame", node_name, item.uri, args)
        _bannerize(args, "count", item.count)
        _bannerize(args, "domain", item.domain)
        _bannerize(args, "maxdomain", item.maxdomain)
        _bannerize(args, "upgraded", item.tiledbsoma_has_upgraded_domain)

    elif isinstance(item, tiledbsoma.SparseNDArray):
        _print_leaf_node_banner("SparseNDArray", node_name, item.uri, args)
        _bannerize(args, "used_shape", item.used_shape())
        _bannerize(args, "shape",      item.shape)
        _bannerize(args, "maxshape",   item.maxshape)
        _bannerize(args, "upgraded",   item.tiledbsoma_has_upgraded_shape)

    elif isinstance(item, tiledbsoma.DenseNDArray):
        _print_leaf_node_banner("DenseNDArray", node_name, item.uri, args)
        _bannerize(args, "shape",      item.shape)
        _bannerize(args, "maxshape",   item.maxshape)
        _bannerize(args, "upgraded",   item.tiledbsoma_has_upgraded_shape)

    return retval


def _leaf_visitor_upgrade(
    item: Any,
    *,
    node_name: str,
    args: SizingArgs,
) -> bool:
    verbose = args["verbose"]
    check_only = args["check_only"]
    retval = True

    if isinstance(item, tiledbsoma.DataFrame):
        count = item.count

        _print_leaf_node_banner("DataFrame", node_name, item.uri, args)
        if check_only:
            print(f"  Dry run for: upgrade_soma_joinid_shape({count})")
            ok, msg = item.upgrade_soma_joinid_shape(count, check_only=True)
            _print_dry_run_result(ok: bool, msg: str)
            retval = retval and ok
        elif not item.tiledbsoma_has_upgraded_domain:
            if verbose:
                print(f"  Applying upgrade_soma_joinid_shape({count})")
            with tiledbsoma.DataFrame.open(item.uri, "w") as writer:
                writer.upgrade_soma_joinid_shape(count)
        else:
            if verbose:
                print("  Already upgraded")

    elif isinstance(item, tiledbsoma.SparseNDArray):
        used_shape = item.used_shape()
        new_shape = tuple(e[1] + 1 for e in used_shape)

        _print_leaf_node_banner("SparseNDArray", node_name, item.uri, args)
        if check_only:
            print(f"  Dry run for: upgrade_shape({new_shape})")
            ok, msg = item.tiledbsoma_upgrade_shape(new_shape, check_only=True)
            _print_dry_run_result(ok: bool, msg: str)
            retval = retval and ok
        elif not item.tiledbsoma_has_upgraded_shape:
            if verbose:
                print(f"  Applying upgrade_shape({new_shape})")
            with tiledbsoma.SparseNDArray.open(item.uri, "w") as writer:
                writer.tiledbsoma_upgrade_shape(new_shape)
        else:
            if verbose:
                print("  Already upgraded")

    elif isinstance(item, tiledbsoma.DenseNDArray):
        _print_leaf_node_banner("DenseNDArray", node_name, item.uri, args)
        if verbose or check_only:
            print(
                "  No action at this time, pending new-shape support for dense arrays in core 2.27"
            )
    return retval


def _leaf_visitor_resize(
    item: Any,
    *,
    node_name: str,
    args: SizingArgs,
) -> bool:
    verbose = args["verbose"]
    check_only = args["check_only"]
    retval = True
    if isinstance(item, tiledbsoma.DataFrame):
        nobs = args["nobs"]
        if nobs is None:
            raise tiledbsoma.SOMAError(
                "experiment resize: internal error: nobs missing"
            )

        _print_leaf_node_banner("DataFrame", node_name, item.uri, args)

        if check_only:
            print(f"  Dry run for: resize_soma_joinid_shape({nobs})")
            ok, msg = item.resize_soma_joinid_shape(nobs, check_only=True)
            _print_dry_run_result(ok: bool, msg: str)
            retval = retval and ok
            XXX FINISH TYPING
        elif not item.tiledbsoma_has_upgraded_domain:
            if verbose:
                print(f"  Applying resize_soma_joinid_shape({nobs})")
            with tiledbsoma.DataFrame.open(item.uri, "w") as writer:
                writer.resize_soma_joinid_shape(nobs)
        else:
            if verbose:
                print("  Already upgraded")

    elif isinstance(item, tiledbsoma.SparseNDArray):

        _print_leaf_node_banner("SparseNDArray", node_name, item.uri, args)

        new_shape = _get_new_ndarray_shape(args, item.shape)

        if check_only:
            print(f"  Dry run for: resize({new_shape})")
            ok, msg = item.resize(new_shape, check_only=True)
            if ok:
                print("  OK")
            else:
                print(f"  {msg}")
                retval = False
        elif not item.tiledbsoma_has_upgraded_shape:
            if verbose:
                print(f"  Applying resize({new_shape})")
            with tiledbsoma.SparseNDArray.open(item.uri, "w") as writer:
                writer.resize(new_shape)
        else:
            if verbose:
                print("  Already upgraded")

    elif isinstance(item, tiledbsoma.DenseNDArray):

        _print_leaf_node_banner("DenseNDArray", node_name, item.uri, args)
        if verbose or check_only:
            print(
                "  No action at this time, pending new-shape support for dense arrays in core 2.27"
            )

    return retval


def _print_leaf_node_banner(
    type_name: str,
    node_name: str,
    uri: str,
    args: SizingArgs,
) -> None:
    if args["verbose"] or args["check_only"]:
        print()
        print(_get_leaf_node_description(type_name, node_name, args))
        print(f"  URI {uri}")


def _get_leaf_node_description(
    type_name: str,
    node_name: str,
    args: SizingArgs,
) -> str:
    # Return things like "ms/RNA/X/data"
    pieces = []
    if args["ms_name"] is not None:
        pieces.append("ms")
        pieces.append(args["ms_name"])
    if args["coll_name"] is not None:
        pieces.append(args["coll_name"])
    pieces.append(node_name)
    return "[" + type_name + "] " + "/".join(pieces)

def _bannerize(
    args: SizingArgs,
    name: str,
    value: Any,
) -> None:

    if args["verbose"] or args["check_only"]:
        padded_name = "%-20s" % name
        print(f"  {padded_name} {value}")

def _print_dry_run_result(ok: bool, msg: str) -> None:
    if ok:
        print("  OK")
    else:
        print(f"  {msg}")
        retval = False

def _get_new_ndarray_shape(
    args: SizingArgs,
    current_shape: Tuple[int, ...],
) -> Tuple[int, ...]:

    nobs = args["nobs"]
    if nobs is None:
        raise tiledbsoma.SOMAError("experiment resize: internal error: nobs missing")

    nvars = args["nvars"]
    if nvars is None:
        raise tiledbsoma.SOMAError("experiment resize: internal error: nvars missing")

    ms_name = args["ms_name"]
    if ms_name is None:
        raise tiledbsoma.SOMAError("experiment resize: internal error: ms_name missing")

    if ms_name not in nvars:
        raise tiledbsoma.SOMAError(
            f"experiment resize: missing measurement name '{ms_name}' in provided nvars"
        )
    nvar = nvars[ms_name]

    coll_name = args["coll_name"]
    if coll_name is None:
        raise tiledbsoma.SOMAError(
            "experiment resize: internal error: coll_name missing"
        )

    if coll_name == "X":
        return (nobs, nvar)
    elif coll_name == "obsm":
        return (nobs, current_shape[1])
    elif coll_name == "obsp":
        return (nobs, nobs)
    elif coll_name == "varm":
        return (nvar, current_shape[1])
    elif coll_name == "varp":
        return (nvar, nvar)
    else:
        raise tiledbsoma.SOMAError(
            f"experiment resize: internal error: unhandled collection {coll_name}"
        )
