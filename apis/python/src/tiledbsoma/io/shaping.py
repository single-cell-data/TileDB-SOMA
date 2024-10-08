# Copyright (c) 2024 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2024 TileDB, Inc.
#
# Licensed under the MIT License.

"""
Provides support for the new shape feature in TileDB-SOMA 1.15, including the
ability to process all dataframes/arrays contained within a TileDB-SOMA
Experiment. Please also see
https://github.com/single-cell-data/TileDB-SOMA/issues/2407.  """

from typing import Any, Dict, Optional

import tiledbsoma


def show_experiment_shapes(
    uri: str,
    *,
    context: Optional[tiledbsoma.SOMATileDBContext] = None,
) -> None:
    """For each dataframe/array contained within the SOMA ``Experiment`` pointed
    to by the given URI, shows the deprecated ``used_shape`` (for N-D arrays) or
    the ``count`` (for dataframes), along with the ``shape`` and ``maxshape``
    (for arrays) or ``domain`` and ``maxdomain`` (for dataframes).
    """
    _treewalk(
        uri,
        visitor=_leaf_visitor_show_shapes,
        verbose=True,
        check_only=True,
        ms_name="N/A",
        nobs=-1,
        nvars={},
        context=context,
    )


def upgrade_experiment_shapes(
    uri: str,
    *,
    verbose: bool = False,
    check_only: bool = False,
    context: Optional[tiledbsoma.SOMATileDBContext] = None,
) -> None:
    """For each dataframe contained within the SOMA ``Experiment`` pointed to by
    the given URI, sets the ``domain`` to match the dataframe's current
    ``count``.  For each N-D array, sets the ``shape`` to match the array's
    current ``used_shape``. If ``verbose`` is set to ``True``, an activity log
    is printed. If ``check_only`` is true, only does a dry run and reports any
    reasons the upgrade would fail.
    """
    _treewalk(
        uri,
        visitor=_leaf_visitor_upgrade,
        verbose=verbose,
        check_only=check_only,
        ms_name="N/A",
        nobs=-1,
        nvars={},
        context=context,
    )
    pass


def resize_experiment(
    uri: str,
    *,
    verbose: bool = False,
    check_only: bool = False,
    nobs: int,
    nvars: Dict[str, int],
    context: Optional[tiledbsoma.SOMATileDBContext] = None,
) -> None:
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
    """
    pass


def _treewalk(
    uri: str,
    *,
    visitor: Any,
    verbose: bool,
    check_only: bool,
    ms_name: str,
    nobs: int,
    nvars: Dict[str, int],
    context: Optional[tiledbsoma.SOMATileDBContext],
) -> None:
    with tiledbsoma.open(uri) as item:
        if isinstance(item, tiledbsoma.Experiment):
            if "obs" in item:
                _treewalk(
                    item["obs"].uri,
                    visitor=visitor,
                    verbose=verbose,
                    check_only=check_only,
                    ms_name=ms_name,
                    nobs=nobs,
                    nvars=nvars,
                    context=context,
                )
            if "ms" in item:
                _treewalk(
                    item["ms"].uri,
                    visitor=visitor,
                    verbose=verbose,
                    check_only=check_only,
                    ms_name=ms_name,
                    nobs=nobs,
                    nvars=nvars,
                    context=context,
                )

        elif isinstance(item, tiledbsoma.Measurement):
            for coll_name in ["X", "obsm", "obsp", "varm", "varp"]:
                if coll_name in item:
                    _treewalk(
                        item[coll_name].uri,
                        visitor=visitor,
                        verbose=verbose,
                        check_only=check_only,
                        ms_name=ms_name,
                        nobs=nobs,
                        nvars=nvars,
                        context=context,
                    )

        elif isinstance(item, tiledbsoma.Collection):
            for key in item:
                _treewalk(
                    item[key].uri,
                    visitor=visitor,
                    verbose=verbose,
                    check_only=check_only,
                    ms_name=ms_name,
                    nobs=nobs,
                    nvars=nvars,
                    context=context,
                )

        else:
            visitor(
                item,
                verbose=verbose,
                check_only=check_only,
                ms_name=ms_name,
                nobs=nobs,
                nvars=nvars,
                context=context,
            )


def _leaf_visitor_show_shapes(
    item: Any,
    *,
    verbose: bool,  # ignored here
    check_only: bool,  # ignored here
    ms_name: str,
    nobs: int,
    nvars: Dict[str, int],
    context: Optional[tiledbsoma.SOMATileDBContext] = None,
) -> None:
    if isinstance(item, tiledbsoma.DataFrame):
        print()
        print("DataFrame")
        print(f"  URI        {item.uri}")
        print(f"  count      {item.count}")
        print(f"  domain     {item.domain}")
        print(f"  maxdomain  {item.maxdomain}")
        print(f"  upgraded   {item.tiledbsoma_has_upgraded_domain}")

    elif isinstance(item, tiledbsoma.SparseNDArray):
        print()
        print("SparseNDArray")
        print(f"  URI        {item.uri}")
        print(f"  used_shape {item.used_shape()}")
        print(f"  shape      {item.shape}")
        print(f"  maxshape   {item.maxshape}")
        print(f"  upgraded   {item.tiledbsoma_has_upgraded_shape}")

    elif isinstance(item, tiledbsoma.DenseNDArray):
        print()
        print("DenseNDArray")
        print(f"  URI        {item.uri}")
        print(f"  shape      {item.shape}")
        print(f"  maxshape   {item.maxshape}")
        print(f"  upgraded   {item.tiledbsoma_has_upgraded_shape}")


def _leaf_visitor_upgrade(
    item: Any,
    *,
    verbose: bool,
    check_only: bool,
    ms_name: str = "",
    nobs: int,
    nvars: Dict[str, int],
    context: Optional[tiledbsoma.SOMATileDBContext] = None,
) -> None:
    if isinstance(item, tiledbsoma.DataFrame):
        count = item.count

        if verbose or check_only:
            print()
            print("DataFrame")
            print(f"  URI        {item.uri}")

        if check_only:
            print(f"  Dry run for: upgrade_soma_joinid_shape({count})")
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

        if verbose or check_only:
            print()
            print("SparseNDArray")
            print(f"  URI        {item.uri}")

        if check_only:
            print(f"  Dry run for: upgrade_shape({new_shape})")
        elif not item.tiledbsoma_has_upgraded_shape:
            if verbose:
                print(f"  Applying upgrade_shape({new_shape})")
            with tiledbsoma.SparseNDArray.open(item.uri, "w") as writer:
                writer.tiledbsoma_upgrade_shape(new_shape)
        else:
            if verbose:
                print("  Already upgraded")

    elif isinstance(item, tiledbsoma.SparseNDArray):

        if verbose or check_only:
            print()
            print("DenseNDArray")
            print(f"  URI        {item.uri}")
            print(
                "  No action at this time, pending new-shape support for dense arrays in core 2.27"
            )


def _leaf_visitor_resize(
    item: Any,
    *,
    verbose: bool,
    check_only: bool,
    ms_name: str = "",
    nobs: int,
    nvars: Dict[str, int],
    context: Optional[tiledbsoma.SOMATileDBContext] = None,
) -> None:
    # XXX NEED MS NAME WIRED THROUGH
    if isinstance(item, tiledbsoma.DataFrame):
        count = item.count

        if verbose or check_only:
            print()
            print("DataFrame")
            print(f"  URI        {item.uri}")

        if check_only:
            print(f"  Dry run for: upgrade_soma_joinid_shape({count})")
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

        if verbose or check_only:
            print()
            print("SparseNDArray")
            print(f"  URI        {item.uri}")

        if check_only:
            print(f"  Dry run for: upgrade_shape({new_shape})")
        elif not item.tiledbsoma_has_upgraded_shape:
            if verbose:
                print(f"  Applying upgrade_shape({new_shape})")
            with tiledbsoma.SparseNDArray.open(item.uri, "w") as writer:
                writer.tiledbsoma_upgrade_shape(new_shape)
        else:
            if verbose:
                print("  Already upgraded")

    elif isinstance(item, tiledbsoma.SparseNDArray):

        if verbose or check_only:
            print()
            print("DenseNDArray")
            print(f"  URI        {item.uri}")
            print(
                "  No action at this time, pending new-shape support for dense arrays in core 2.27"
            )
