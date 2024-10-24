# Copyright (c) 2024 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2024 TileDB, Inc.
#
# Licensed under the MIT License.

"""
Provides support for the new shape feature in TileDB-SOMA 1.15, including the
ability to process all dataframes/arrays contained within a TileDB-SOMA
Experiment. Please also see
https://github.com/single-cell-data/TileDB-SOMA/issues/2407.  """

import io
import sys
from typing import Any, Dict, Optional, Tuple, TypedDict, Union, cast

import tiledbsoma

Printable = Union[io.TextIOWrapper, io.StringIO]


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
    output_handle: Printable


def show_experiment_shapes(
    uri: str,
    *,
    context: Optional[tiledbsoma.SOMATileDBContext] = None,
    output_handle: Printable = cast(Printable, sys.stdout),
) -> bool:
    """For each dataframe/array contained within the SOMA ``Experiment`` pointed
    to by the given URI, shows the deprecated ``used_shape`` (for N-D arrays) or
    the ``count`` (for dataframes), along with the ``shape`` and ``maxshape``
    (for arrays) or ``domain`` and ``maxdomain`` (for dataframes).

    Args:
        uri: The URI of a SOMA :class:`Experiment``.
        context: Optional :class:`SOMATileDBContext``.
    """
    args: SizingArgs = dict(
        nobs=None,
        nvars=None,
        ms_name=None,
        coll_name=None,
        verbose=True,
        check_only=True,
        context=context,
        output_handle=output_handle,
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
    output_handle: Printable = cast(Printable, sys.stdout),
) -> bool:
    """For each dataframe contained within the SOMA ``Experiment`` pointed to by
    the given URI, sets the ``domain`` to match the dataframe's current
    ``count``.  For each N-D array, sets the ``shape`` to match the array's
    current ``used_shape``. If ``verbose`` is set to ``True``, an activity log
    is printed. If ``check_only`` is true, only does a dry run and reports any
    reasons the upgrade would fail.

    Args:
        uri: The URI of a SOMA :class:`Experiment``.
        verbose: If ``True``, produce per-array output as the upgrade runs.
        check_only: If ``True``,  don't apply the upgrades, but show what would
            be attempted, and show why each one would fail.
        context: Optional :class:`SOMATileDBContext``.
    """
    args: SizingArgs = dict(
        nobs=None,
        nvars=None,
        ms_name=None,
        coll_name=None,
        verbose=verbose,
        check_only=check_only,
        context=context,
        output_handle=output_handle,
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
    output_handle: Printable = cast(Printable, sys.stdout),
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

    Args:
        uri: The URI of a SOMA :class:`Experiment``.
        nobs: The desired new shape of the experiment's ``obs`` dataframe.
        nvars: The desired new shapes of the experiment's ``var`` dataframes.
            This should be a dict from measurement name to shape, e.g.
            ``{"RNA": 10000, "raw": 20000}``.
        verbose: If ``True``, produce per-array output as the upgrade runs.
        check_only: If ``True``,  don't apply the upgrades, but show what would
            be attempted, and show why each one would fail.
        context: Optional :class:`SOMATileDBContext``.
    """
    args: SizingArgs = dict(
        nobs=nobs,
        nvars=nvars,
        ms_name=None,
        coll_name=None,
        verbose=verbose,
        check_only=check_only,
        context=context,
        output_handle=output_handle,
    )

    # Extra user-provided keys not relevant to the experiment are ignored.  This
    # is important for the case when a new measurement, which is registered from
    # AnnData/H5AD inputs, is registered and is about to be created but does not
    # exist just yet in the experiment storage.
    #
    # If the user hasn't provided a key -- e.g. a from-anndata-append-with-resize
    # on one measurement while the experiment's other measurements aren't being
    # updated -- then we need to find those other measurements' var-shapes.
    with tiledbsoma.Experiment.open(uri) as exp:
        for ms_key in exp.ms.keys():
            if ms_key not in nvars.keys():
                nvars[ms_key] = exp.ms[ms_key].var._maybe_soma_joinid_shape or 1

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

                    args["coll_name"] = coll_name

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
                    args["ms_name"] = key

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
        _bannerize(args, "shape", item.shape)
        _bannerize(args, "maxshape", item.maxshape)
        _bannerize(args, "upgraded", item.tiledbsoma_has_upgraded_shape)

    elif isinstance(item, tiledbsoma.DenseNDArray):
        _print_leaf_node_banner("DenseNDArray", node_name, item.uri, args)
        _bannerize(args, "shape", item.shape)
        _bannerize(args, "maxshape", item.maxshape)
        _bannerize(args, "upgraded", item.tiledbsoma_has_upgraded_shape)

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
            print(
                f"  Dry run for: tiledbsoma_upgrade_soma_joinid_shape({count})",
                file=args["output_handle"],
            )
            ok, msg = item.tiledbsoma_upgrade_soma_joinid_shape(count, check_only=True)
            _print_dry_run_result(ok, msg, args)
            retval = retval and ok
        elif not item.tiledbsoma_has_upgraded_domain:
            if verbose:
                print(
                    f"  Applying tiledbsoma_upgrade_soma_joinid_shape({count})",
                    file=args["output_handle"],
                )
            with tiledbsoma.DataFrame.open(item.uri, "w") as writer:
                writer.tiledbsoma_upgrade_soma_joinid_shape(count)
        else:
            if verbose:
                print("  Already upgraded", file=args["output_handle"])

    elif isinstance(item, tiledbsoma.SparseNDArray):
        used_shape = item.used_shape()
        new_shape = tuple(e[1] + 1 for e in used_shape)

        _print_leaf_node_banner("SparseNDArray", node_name, item.uri, args)
        if check_only:
            print(
                f"  Dry run for: tiledbsoma_upgrade_shape({new_shape})",
                file=args["output_handle"],
            )
            ok, msg = item.tiledbsoma_upgrade_shape(new_shape, check_only=True)
            _print_dry_run_result(ok, msg, args)
            retval = retval and ok
        elif not item.tiledbsoma_has_upgraded_shape:
            if verbose:
                print(
                    f"  Applying tiledbsoma_upgrade_shape({new_shape})",
                    file=args["output_handle"],
                )
            with tiledbsoma.SparseNDArray.open(item.uri, "w") as writer:
                writer.tiledbsoma_upgrade_shape(new_shape)
        else:
            if verbose:
                print("  Already upgraded", file=args["output_handle"])

    elif isinstance(item, tiledbsoma.DenseNDArray):
        _print_leaf_node_banner("DenseNDArray", node_name, item.uri, args)
        if verbose or check_only:
            print(
                "  No action at this time, pending new-shape support for dense arrays in core 2.27",
                file=args["output_handle"],
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

        if node_name == "obs":
            new_soma_joinid_shape = args["nobs"]
            if new_soma_joinid_shape is None:
                raise tiledbsoma.SOMAError(
                    "experiment resize: internal error: nobs missing"
                )

        elif node_name == "var":
            new_soma_joinid_shape = _get_new_var_shape(args)

        else:
            raise tiledbsoma.SOMAError(
                "experiment resize: internal error: dataframe node name '{node_name}'"
            )

        _print_leaf_node_banner("DataFrame", node_name, item.uri, args)

        if check_only:
            print(
                f"  Dry run for: tiledbsoma_resize_soma_joinid_shape({new_soma_joinid_shape})",
                file=args["output_handle"],
            )
            ok, msg = item.tiledbsoma_resize_soma_joinid_shape(
                new_soma_joinid_shape, check_only=True
            )
            _print_dry_run_result(ok, msg, args)
            retval = retval and ok
        else:
            if verbose:
                print(
                    f"  Applying tiledbsoma_resize_soma_joinid_shape({new_soma_joinid_shape})",
                    file=args["output_handle"],
                )
            with tiledbsoma.DataFrame.open(item.uri, "w") as writer:
                writer.tiledbsoma_resize_soma_joinid_shape(new_soma_joinid_shape)

    elif isinstance(item, tiledbsoma.SparseNDArray):

        _print_leaf_node_banner("SparseNDArray", node_name, item.uri, args)

        new_shape = _get_new_ndarray_shape(args, item.shape)

        if check_only:
            print(f"  Dry run for: resize({new_shape})", file=args["output_handle"])
            ok, msg = item.resize(new_shape, check_only=True)
            _print_dry_run_result(ok, msg, args)
            retval = retval and ok
        else:
            if verbose:
                print(f"  Applying resize({new_shape})", file=args["output_handle"])
            with tiledbsoma.SparseNDArray.open(item.uri, "w") as writer:
                writer.resize(new_shape)

    elif isinstance(item, tiledbsoma.DenseNDArray):

        _print_leaf_node_banner("DenseNDArray", node_name, item.uri, args)
        if verbose or check_only:
            print(
                "  No action at this time, pending new-shape support for dense arrays in core 2.27",
                file=args["output_handle"],
            )

    return retval


def _print_leaf_node_banner(
    type_name: str,
    node_name: str,
    uri: str,
    args: SizingArgs,
) -> None:
    if args["verbose"] or args["check_only"]:
        print("", file=args["output_handle"])
        print(
            _get_leaf_node_description(type_name, node_name, args),
            file=args["output_handle"],
        )
        print(f"  URI {uri}", file=args["output_handle"])


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
    return f"[{type_name}] {'/'.join(pieces)} "


def _bannerize(
    args: SizingArgs,
    name: str,
    value: Any,
) -> None:
    if args["verbose"] or args["check_only"]:
        print(f"  {name:<20} {value}", file=args["output_handle"])


def _print_dry_run_result(ok: bool, msg: str, args: SizingArgs) -> None:
    if ok:
        print("  OK", file=args["output_handle"])
    else:
        print(f"  Not OK: {msg}", file=args["output_handle"])


def _get_new_var_shape(
    args: SizingArgs,
) -> int:
    """Maps from experiment-level nobs and per-measurement nvar to the correct resizes
    for var dataframes within a SOMA experiment.
    """
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
    return nvars[ms_name]


def _get_new_ndarray_shape(
    args: SizingArgs,
    current_shape: Tuple[int, ...],
) -> Tuple[int, ...]:
    """Maps from experiment-level nobs and per-measurement nvar to the correct resizes
    for various kinds of n-d array within a SOMA Experiment.
    """
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

    coll_dict = {
        "X": (nobs, nvar),
        "obsm": (nobs, current_shape[1]),
        "obsp": (nobs, nobs),
        "varm": (nvar, current_shape[1]),
        "varp": (nvar, nvar),
    }

    try:
        return coll_dict[coll_name]
    except KeyError:
        raise tiledbsoma.SOMAError(
            f"experiment resize: internal error: unhandled collection {coll_name}"
        )
