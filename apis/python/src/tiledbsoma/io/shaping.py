# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""Provides support for the new shape feature in TileDB-SOMA 1.15, including the
ability to process all dataframes/arrays contained within a TileDB-SOMA
Experiment. Please also see
https://github.com/single-cell-data/TileDB-SOMA/issues/2407.
"""

from __future__ import annotations

import io
import sys
from typing import Any, Callable, TypeVar, Union, cast

import tiledbsoma

from .._soma_object import SOMAObject

Printable = Union[io.TextIOWrapper, io.StringIO]

_SOMAObjectType = TypeVar("_SOMAObjectType", bound=SOMAObject)  # type: ignore[type-arg]


def get_experiment_shapes(
    uri: str,
    *,
    context: tiledbsoma.SOMATileDBContext | None = None,
) -> dict[str, Any]:
    """Returns the current shapes of the elements in the ``Experiment``.

    Similar to ``show_experiment_shapes`` but returns a nested Python dict.

    Returns the ``non_empty_domain``, ``domain``, and ``maxdomain`` of dataframes and
    the ``non_empty_domain``, ``shape``, and ``maxshape`` of arrays as a dictionary. This
    method is applied to the following elements inside the SOMA ``Experiment``:

    The shapes of the following elements are output:

      * the ``obs`` dataframe in the experiment,

    and for each measurement:

      * the ``var`` dataframe,
      * all ``X`` arrays,
      * all ``obsm`` arrays,
      * all ``varm`` arrays,
      * all ``obsp`` arrays,
      * all ``varm`` arrays,
      * all ``varp`` arrays.

    Example::

        >>> tiledbsoma.io.get_experiment_shapes('pbmc3k_unprocessed')

        {
            "obs": {
                "uri": "exp/obs",
                "type": "DataFrame",
                "count": 2700,
                "non_empty_domain": ((0, 2699),),
                "domain": ((0, 2699),),
                "maxdomain": ((0, 2147483646),),
                "upgraded": False,
            },
            "ms": {
                "RNA": {
                    "var": {
                        "uri": "exp/ms/RNA/var",
                        "type": "DataFrame",
                        "count": 13714,
                        "non_empty_domain": ((0, 13713),),
                        "domain": ((0, 13713),),
                        "maxdomain": ((0, 2147483646),),
                        "upgraded": False,
                    },
                    "X": {
                        "data": {
                            "uri": "exp/ms/RNA/X/data",
                            "type": "SparseNDArray",
                            "non_empty_domain": ((0, 2699), (0, 13713)),
                            "shape": (2700, 13714),
                            "maxshape": (2147483646, 2147483646),
                            "upgraded": False,
                        }
                    },
                }
            },
        }


    Args:
        uri: The URI of a SOMA :class:`Experiment`.
        context: Optional :class:`SOMATileDBContext`.
        output_handle: The handle to print the output to.

    Returns: a nested Python dict.
    """
    # Developer note: it would be a well-intentioned mistake to refactor
    # show_experiment_shapes to call get_experiment_shapes and then print the
    # result. The crucial issue is that for complex experiments on remote object
    # stores, show_experiment_shapes communicates liveness UX to the user by
    # printing things as it finds them, rather than giving the user a wait
    # followed by a final burst of output.

    retval = _treewalk(
        uri,
        leaf_visitor=_leaf_visitor_get_shapes,
        nobs=None,
        nvars=None,
        ms_name=None,
        coll_name=None,
        verbose=True,
        check_only=True,
        context=context,
        output_handle=None,
    )
    return retval


def show_experiment_shapes(
    uri: str,
    *,
    context: tiledbsoma.SOMATileDBContext | None = None,
    output_handle: Printable = cast(Printable, sys.stdout),
) -> bool:
    """Outputs the current shapes of the elements in the ``Experiment``.

    Similar to ``get_experiment_shapes`` but prints results as they are
    traversed.

    Outputs the ``non_empty_domain``, ``domain``, and ``maxdomain`` of dataframes and
    the ``non_empty_domain``, ``shape``, and ``maxshape`` of arrays. This method is
    applied to the following elements inside the SOMA ``Experiment``:

    The shapes of the following elements are output:

      * the ``obs`` dataframe in the experiment,

    and for each measurement:

      * the ``var`` dataframe,
      * all ``X`` arrays,
      * all ``obsm`` arrays,
      * all ``varm`` arrays,
      * all ``obsp`` arrays,
      * all ``varm`` arrays,
      * all ``varp`` arrays.


    Example::

        >>> tiledbsoma.io.show_experiment_shapes('pbmc3k_unprocessed')
        [DataFrame] obs
          URI file:///data/pbmc3k_unprocessed/obs
          non_empty_domain     ((0, 2699),)
          domain               ((0, 2699),)
          maxdomain            ((0, 9223372036854773758),)
          upgraded             True
        [DataFrame] ms/RNA/var
          URI file:///data/pbmc3k_unprocessed/ms/RNA/var
          non_empty_domain     ((0, 13713),)
          domain               ((0, 13713),)
          maxdomain            ((0, 9223372036854773758),)
          upgraded             True
        [SparseNDArray] ms/RNA/X/data
          URI file:///data/pbmc3k_unprocessed/ms/RNA/X/data
          shape                (2700, 13714)
          maxshape             (9223372036854773759, 9223372036854773759)
          upgraded             True

    Args:
        uri: The URI of a SOMA :class:`Experiment`.
        context: Optional :class:`SOMATileDBContext`.
        output_handle: The handle to print the output to.

    Returns:
        ``True`` if outputting the shape works for elements. ``False`` if any element
        fails to successfully output its shape.
    """
    retval = _treewalk(
        uri,
        leaf_visitor=_leaf_visitor_show_shapes,
        nobs=None,
        nvars=None,
        ms_name=None,
        coll_name=None,
        verbose=True,
        check_only=True,
        output_handle=output_handle,
        context=context,
    )
    return _check_statuses(retval)


def upgrade_experiment_shapes(
    uri: str,
    *,
    verbose: bool = False,
    check_only: bool = False,
    context: tiledbsoma.SOMATileDBContext | None = None,
    output_handle: Printable = cast(Printable, sys.stdout),
) -> bool:
    """Upgrade the elements inside a SOMA ``Experiment`` to use the ``shape`` feature
    introduced in TileDB-SOMA 1.15.

    A new shape feature was introduced in TileDB-SOMA in release 1.15. This updates
    the elements in TileDB-SOMA to use the new feature. It makes an experiment created
    before TileDB-SOMA 1.15 look like an experiment created by TileDB-SOMA 1.15 or
    later. You can use ``tiledbsoma.io.show_experiment_shapes`` before and after to see
    the difference.

    For each dataframe and N-D array that is being upgraded, if the dataframe does
    not currently support the new shape feature, upgrades to add the feature and sets
    the domain to the dataframe's current ``non_empty_domain``. If the new ``shape``
    feature is already enable, nothing is changed.

    The following elements are updated:

      * the ``obs`` dataframe in the experiment,

    for each measurement:

      * the ``var`` dataframe,
      * all ``X`` arrays,
      * all ``obsm`` arrays,
      * all ``varm`` arrays,
      * all ``obsp`` arrays,
      * all ``varm`` arrays,
      * all ``varp`` arrays.

    Example::

        >>> tiledbsoma.io.upgrade_experiment_shapes('pbmc3k_unprocessed_old', check_only=True)
        [DataFrame] obs
          URI file:///data/pbmc3k_unprocessed_old/obs
          Dry run for: tiledbsoma_upgrade_soma_joinid_shape(2700)
          OK
        [DataFrame] ms/RNA/var
          URI file:///data/pbmc3k_unprocessed_old/ms/RNA/var
          Dry run for: tiledbsoma_upgrade_soma_joinid_shape(13714)
          OK
        [SparseNDArray] ms/RNA/X/data
          URI file:///data/pbmc3k_unprocessed_old/ms/RNA/X/data
          Dry run for: tiledbsoma_upgrade_shape((2700, 13714))
          OK

    Args:
        uri: The URI of a SOMA :class:`Experiment`.
        verbose: If ``True``, produce per-array output as the upgrade runs.
        check_only: If ``True``,  don't apply the upgrades, but show what would
            be attempted, and show why each one would fail.
        context: Optional :class:`SOMATileDBContext`.

    Returns:
        ``True`` if all upgrade operations succeed. ``False`` if any upgrade
        operation fails.
    """
    # For resize, nobs and nvars are from the user. But for upgrade,
    # they're from the experiment as-is.
    nobs = None
    nvars = {}
    with tiledbsoma.Experiment.open(uri, context=context) as exp:
        if "obs" in exp:
            nobs = exp.obs.count
        if "ms" in exp:
            for name, measurement in exp.ms.items():
                if "var" in measurement:
                    nvars[name] = measurement.var.count

    retval = _treewalk(
        uri,
        leaf_visitor=_leaf_visitor_upgrade,
        nobs=nobs,
        nvars=nvars,
        ms_name=None,
        coll_name=None,
        verbose=verbose,
        check_only=check_only,
        output_handle=output_handle,
        context=context,
    )
    return _check_statuses(retval)


def resize_experiment(
    uri: str,
    *,
    nobs: int,
    nvars: dict[str, int],
    verbose: bool = False,
    check_only: bool = False,
    context: tiledbsoma.SOMATileDBContext | None = None,
    output_handle: Printable = cast(Printable, sys.stdout),
) -> bool:
    """Resize the elements in the SOMA ``Experiment`` to fit the requested number
    of observations and variables.

    A dataframe will be resized if the ``soma_joinid`` column is an index column
    and the current domain of the ``soma_joinid`` column is smaller than requested
    size.  If the new domain of the ``soma_joinid`` column does not fit inside the
    ``maxshape``, the resize fails.

    An N-D array will be resized if either dimension of the new ``shape`` is larger
    than the current shape. If either dimension of the new ``shape`` is larger
    than ``maxshape``, the resize fails.

    The following base elements are resized:

      * ``obs`` dataframe: ``soma_joinid`` resized to fit ``nobs``.

    For each ``measurement_name`` in the experiment the elements are resized as follows
    where ``nvar[measurement_name]`` is the current size if ``measurement_name`` or
    no value is provided by the user:

      * ``var`` dataframe: ``soma_joinid`` resized to fit ``nvar[measurement_name]``.
      * ``X`` arrays: ``shape`` at least (``nobs``, ``nvar[measurement_name]``).
      * ``obsm`` arrays: ``shape`` at least (``nobs``, current ``soma_dim_1``).
      * ``varm`` arrays: ``shape`` at least (``nvar``, current ``soma_dim_1``).
      * ``obsp`` arrays: ``shape`` at least (``nobs``, ``nobs``).
      * ``varm`` arrays: ``shape`` at least (``nvar``, existing ``soma_dim_1``).
      * ``varp`` arrays: ``shape`` at least ( ``nvar``, ``nvar``).

    Example::

        >>> tiledbsoma.io.resize_experiment(
            'pbmc3k_unprocessed',
            nobs=5600,
            nvars={"data":3204, "raw": 13714},
            check_only=True,
        )
        [DataFrame] obs
          URI file:///data/pbmc3k_unprocessed/obs
          Dry run for: tiledbsoma_resize_soma_joinid_shape(5600)
          OK
        [DataFrame] ms/RNA/var
          URI file:///data/pbmc3k_unprocessed/ms/RNA/var
          Dry run for: tiledbsoma_resize_soma_joinid_shape(13714)
          OK
        [SparseNDArray] ms/RNA/X/data
          URI file:///data/pbmc3k_unprocessed/ms/RNA/X/data
          Dry run for: resize((5600, 13714))
          OK


    Args:
        uri: The URI of a SOMA :class:`Experiment`.
        nobs: The desired new shape of the experiment's ``obs`` dataframe.
        nvars: The desired new shapes of the experiment's ``var`` dataframes.
            This should be a dict from measurement name to shape, e.g.
            ``{"RNA": 10000, "raw": 20000}``.
        verbose: If ``True``, produce per-array output as the upgrade runs.
        check_only: If ``True``,  don't apply the upgrades, but show what would
            be attempted, and show why each one would fail.
        context: Optional :class:`SOMATileDBContext`.

    Returns:
        ``True`` if all resize operations succeed. ``False`` if any resize operation
        fails.
    """
    # Extra user-provided keys not relevant to the experiment are ignored.  This
    # is important for the case when a new measurement, which is registered from
    # AnnData/H5AD inputs, is registered and is about to be created but does not
    # exist just yet in the experiment storage.
    #
    # If the user hasn't provided a key -- e.g. a from-anndata-append-with-resize
    # on one measurement while the experiment's other measurements aren't being
    # updated -- then we need to find those other measurements' var-shapes.
    with tiledbsoma.Experiment.open(uri, context=context) as exp:
        for ms_key in exp.ms.keys():
            if ms_key not in nvars.keys():
                nvars[ms_key] = exp.ms[ms_key].var._maybe_soma_joinid_shape or 1

    retval = _treewalk(
        uri,
        leaf_visitor=_leaf_visitor_resize,
        nobs=nobs,
        nvars=nvars,
        ms_name=None,
        coll_name=None,
        verbose=verbose,
        check_only=check_only,
        output_handle=output_handle,
        context=context,
    )
    return _check_statuses(retval)


def _treewalk(
    uri: str,
    *,
    node_name: str | None = None,
    leaf_visitor: Callable[..., dict[str, Any]],
    **kwargs: Any,
) -> dict[str, Any]:
    """Apply visitor function to the ``Experiment`` elements.

    The following elements are visited:

      * the ``obs`` dataframe in the experiment,

    for each measurement:

      * the ``var`` dataframe,
      * all ``X`` arrays,
      * all ``obsm`` arrays,
      * all ``varm`` arrays,
      * all ``obsp`` arrays,
      * all ``varm`` arrays,
      * all ``varp`` arrays.

    Args:
        uri: URI of the element to visit and visit the children of.
        node_name: Name of the element to visit and visit the children of.
        leaf_visitor: A callback function that returns a dict for a given
            leaf node (DataFrame or ND Array), along with doing any
            side effects such as printing to the output handle.

        Nominal kwargs:

        nobs: int | None: Gotten from Experiment obs and passed into Experiment ms.
        nvars: Dict[str, int] | None: Gotten from Measurement var and passed into
            Measurement collections such as X, obsm, etc.
        ms_name: str | None: Measurement name for elements of Experiment ms.
        coll_name: str | None: "X", "obsm", etc. passed to members of those collections.
        verbose: bool: Tracks user-level flags.
        check_only: bool: Tracks user-level flags.
        output_handle: Printable | None: Nominally sys.stdout.
        context: TileDB context to pass to the visitor.
    """

    def _recurse(
        parent: (
            tiledbsoma.Experiment
            | tiledbsoma.Measurement
            | tiledbsoma.Collection[_SOMAObjectType]
        ),
        node_name: str | None,
        **kwargs: Any,
    ) -> dict[str, Any]:
        """Applies ``_treewalk`` to the requested child element.

        Args:
            parent: The parent item to get child from.
            name: Name of the child element to attempt to visit.
            kwargs: as in ``_treewalk``.

        Returns:
            A nested map with results from all leaf nodes.
        """
        if node_name is not None and node_name in parent:
            return _treewalk(
                uri=parent[node_name].uri,
                node_name=node_name,
                leaf_visitor=leaf_visitor,
                **kwargs,
            )
        return {}

    # Don't abort the treewalk, except in case of exceptions of course.
    # When users want to see if an experiment is upgradeable, they need
    # to be shown *all* the component arrays which have issues, not just
    # the first such.
    with tiledbsoma.open(uri, context=kwargs["context"]) as item:

        # Non-terminal
        if isinstance(item, tiledbsoma.Experiment):
            retval = {}
            if "obs" in item:
                retval["obs"] = _recurse(item, node_name="obs", **kwargs)
            if "ms" in item:
                retval["ms"] = _recurse(item, node_name="ms", **kwargs)
            return retval

        # Non-terminal
        if isinstance(item, tiledbsoma.Measurement):
            retval = {}
            if "var" in item:
                retval["var"] = _recurse(item, node_name="var", **kwargs)
            for iter_coll_name in ["X", "obsm", "obsp", "varm", "varp"]:
                kwargs["coll_name"] = iter_coll_name
                children_dict = _recurse(item, node_name=iter_coll_name, **kwargs)
                # If the experiment has no varp (etc.) at all, don't include an entry for it.
                if children_dict:
                    retval[iter_coll_name] = children_dict
            return retval

        # Non-terminal
        if isinstance(item, tiledbsoma.Collection):
            retval = {}
            for key in item:
                if node_name == "ms":
                    kwargs["ms_name"] = key
                retval[key] = _recurse(item, node_name=key, **kwargs)
            return retval

        # Terminal
        return leaf_visitor(item, node_name=node_name, **kwargs)


def _leaf_visitor_show_shapes(
    item: Any,
    *,
    node_name: str,
    nobs: int | None,
    nvars: dict[str, int] | None,
    ms_name: str | None,
    coll_name: str | None,
    verbose: bool,
    check_only: bool,
    output_handle: Printable | None,
    context: tiledbsoma.SOMATileDBContext | None,
) -> dict[str, Any]:
    retval = {"status": True}
    if isinstance(item, tiledbsoma.DataFrame):
        _print_leaf_node_banner(
            uri=item.uri,
            type_name="DataFrame",
            node_name=node_name,
            ms_name=ms_name,
            coll_name=coll_name,
            verbose=verbose,
            check_only=check_only,
            output_handle=output_handle,
        )
        _bannerize(
            "count",
            item.count,
            verbose=verbose,
            check_only=check_only,
            output_handle=output_handle,
        )
        _bannerize(
            "non_empty_domain",
            item.non_empty_domain(),
            verbose=verbose,
            check_only=check_only,
            output_handle=output_handle,
        )
        _bannerize(
            "domain",
            item.domain,
            verbose=verbose,
            check_only=check_only,
            output_handle=output_handle,
        )
        _bannerize(
            "maxdomain",
            item.maxdomain,
            verbose=verbose,
            check_only=check_only,
            output_handle=output_handle,
        )
        _bannerize(
            "upgraded",
            item.tiledbsoma_has_upgraded_domain,
            verbose=verbose,
            check_only=check_only,
            output_handle=output_handle,
        )

    elif isinstance(item, tiledbsoma.SparseNDArray):
        _print_leaf_node_banner(
            uri=item.uri,
            type_name="SparseNDArray",
            node_name=node_name,
            ms_name=ms_name,
            coll_name=coll_name,
            verbose=verbose,
            check_only=check_only,
            output_handle=output_handle,
        )
        _bannerize(
            "non_empty_domain",
            item.non_empty_domain(),
            verbose=verbose,
            check_only=check_only,
            output_handle=output_handle,
        )
        _bannerize(
            "shape",
            item.shape,
            verbose=verbose,
            check_only=check_only,
            output_handle=output_handle,
        )
        _bannerize(
            "maxshape",
            item.maxshape,
            verbose=verbose,
            check_only=check_only,
            output_handle=output_handle,
        )
        _bannerize(
            "upgraded",
            item.tiledbsoma_has_upgraded_shape,
            verbose=verbose,
            check_only=check_only,
            output_handle=output_handle,
        )

    elif isinstance(item, tiledbsoma.DenseNDArray):
        _print_leaf_node_banner(
            uri=item.uri,
            type_name="DenseNDArray",
            node_name=node_name,
            ms_name=ms_name,
            coll_name=coll_name,
            verbose=verbose,
            check_only=check_only,
            output_handle=output_handle,
        )
        _bannerize(
            "shape",
            item.shape,
            verbose=verbose,
            check_only=check_only,
            output_handle=output_handle,
        )
        _bannerize(
            "maxshape",
            item.maxshape,
            verbose=verbose,
            check_only=check_only,
            output_handle=output_handle,
        )
        _bannerize(
            "upgraded",
            item.tiledbsoma_has_upgraded_shape,
            verbose=verbose,
            check_only=check_only,
            output_handle=output_handle,
        )

    return retval


def _leaf_visitor_upgrade(
    item: Any,
    *,
    node_name: str,
    nobs: int | None,
    nvars: dict[str, int] | None,
    ms_name: str | None,
    coll_name: str | None,
    verbose: bool,
    check_only: bool,
    output_handle: Printable | None,
    context: tiledbsoma.SOMATileDBContext | None,
) -> dict[str, Any]:
    retval = {"status": True}

    if isinstance(item, tiledbsoma.DataFrame):
        if item.index_column_names == ("soma_joinid",):
            count = item.non_empty_domain()[0][1] + 1
        else:
            count = item.count

        _print_leaf_node_banner(
            uri=item.uri,
            type_name="DataFrame",
            node_name=node_name,
            ms_name=ms_name,
            coll_name=coll_name,
            verbose=verbose,
            check_only=check_only,
            output_handle=output_handle,
        )
        if check_only:
            print(
                f"  Dry run for: tiledbsoma_upgrade_soma_joinid_shape({count})",
                file=output_handle,
            )
            ok, msg = item.tiledbsoma_upgrade_soma_joinid_shape(count, check_only=True)
            _print_dry_run_result(
                ok=ok,
                msg=msg,
                output_handle=output_handle,
            )
            retval["status"] = retval["status"] and ok
        elif not item.tiledbsoma_has_upgraded_domain:
            if verbose:
                print(
                    f"  Applying tiledbsoma_upgrade_soma_joinid_shape({count})",
                    file=output_handle,
                )
            with tiledbsoma.DataFrame.open(item.uri, "w", context=context) as writer:
                writer.tiledbsoma_upgrade_soma_joinid_shape(count)
        else:
            if verbose:
                print("  Already upgraded", file=output_handle)

    elif isinstance(item, tiledbsoma.SparseNDArray):

        old_bounds = item.non_empty_domain()
        # Make a tuple of hi+1
        counts_from_old_bounds = tuple(e[1] + 1 for e in old_bounds)
        # Get the right thing to do for X, obsm, obsp, varm, or varp.  Note that
        # counts_from_old_bounds can be the wrong thing in case a given sparse array
        # has _no_ occupied cells in the last one or more rows.  E.g. if nobs is
        # 1000 and nvars["RNA"] is 200 then each X needs to have shape like 1000
        # x 200 but its old bounds might only be say 998 x 200 if there are no
        # X-counts at the end of X.
        new_shape = _get_new_ndarray_shape(
            current_shape=counts_from_old_bounds,
            nobs=nobs,
            nvars=nvars,
            ms_name=ms_name,
            coll_name=coll_name,
        )

        _print_leaf_node_banner(
            uri=item.uri,
            type_name="SparseNDArray",
            node_name=node_name,
            ms_name=ms_name,
            coll_name=coll_name,
            verbose=verbose,
            check_only=check_only,
            output_handle=output_handle,
        )
        if check_only:
            print(
                f"  Dry run for: tiledbsoma_upgrade_shape({new_shape})",
                file=output_handle,
            )
            ok, msg = item.tiledbsoma_upgrade_shape(new_shape, check_only=True)
            _print_dry_run_result(
                ok=ok,
                msg=msg,
                output_handle=output_handle,
            )
            retval["status"] = retval["status"] and ok
        elif not item.tiledbsoma_has_upgraded_shape:
            if verbose:
                print(
                    f"  Applying tiledbsoma_upgrade_shape({new_shape})",
                    file=output_handle,
                )
            with tiledbsoma.SparseNDArray.open(
                item.uri, "w", context=context
            ) as writer:
                writer.tiledbsoma_upgrade_shape(new_shape)
        else:
            if verbose:
                print("  Already upgraded", file=output_handle)

    elif isinstance(item, tiledbsoma.DenseNDArray):
        _print_leaf_node_banner(
            uri=item.uri,
            type_name="DenseNDArray",
            node_name=node_name,
            ms_name=ms_name,
            coll_name=coll_name,
            verbose=verbose,
            check_only=check_only,
            output_handle=output_handle,
        )
        if verbose or check_only:
            print(
                "  No action at this time, pending new-shape support for dense arrays in core 2.27",
                file=output_handle,
            )
    return retval


def _leaf_visitor_resize(
    item: Any,
    *,
    node_name: str,
    nobs: int | None,
    nvars: dict[str, int] | None,
    ms_name: str | None,
    coll_name: str | None,
    verbose: bool,
    check_only: bool,
    context: tiledbsoma.SOMATileDBContext | None,
    output_handle: Printable | None,
) -> dict[str, Any]:
    retval = {"status": True}
    if isinstance(item, tiledbsoma.DataFrame):

        if node_name == "obs":
            new_soma_joinid_shape = nobs
            if new_soma_joinid_shape is None:
                raise tiledbsoma.SOMAError(
                    "experiment resize: internal error: nobs missing"
                )

        elif node_name == "var":
            new_soma_joinid_shape = _get_new_var_shape(nvars=nvars, ms_name=ms_name)

        else:
            raise tiledbsoma.SOMAError(
                "experiment resize: internal error: dataframe node name '{node_name}'"
            )

        _print_leaf_node_banner(
            uri=item.uri,
            type_name="DataFrame",
            node_name=node_name,
            ms_name=ms_name,
            coll_name=coll_name,
            verbose=verbose,
            check_only=check_only,
            output_handle=output_handle,
        )

        if check_only:
            print(
                f"  Dry run for: tiledbsoma_resize_soma_joinid_shape({new_soma_joinid_shape})",
                file=output_handle,
            )
            ok, msg = item.tiledbsoma_resize_soma_joinid_shape(
                new_soma_joinid_shape, check_only=True
            )
            _print_dry_run_result(
                ok=ok,
                msg=msg,
                output_handle=output_handle,
            )
            retval["status"] = retval["status"] and ok
        else:
            if verbose:
                print(
                    f"  Applying tiledbsoma_resize_soma_joinid_shape({new_soma_joinid_shape})",
                    file=output_handle,
                )
            with tiledbsoma.DataFrame.open(item.uri, "w", context=context) as writer:
                writer.tiledbsoma_resize_soma_joinid_shape(new_soma_joinid_shape)

    elif isinstance(item, tiledbsoma.SparseNDArray):

        _print_leaf_node_banner(
            uri=item.uri,
            type_name="SparseNDArray",
            node_name=node_name,
            ms_name=ms_name,
            coll_name=coll_name,
            verbose=verbose,
            check_only=check_only,
            output_handle=output_handle,
        )

        new_shape = _get_new_ndarray_shape(
            current_shape=item.shape,
            nobs=nobs,
            nvars=nvars,
            ms_name=ms_name,
            coll_name=coll_name,
        )

        if check_only:
            print(f"  Dry run for: resize({new_shape})", file=output_handle)
            ok, msg = item.resize(new_shape, check_only=True)
            _print_dry_run_result(
                ok=ok,
                msg=msg,
                output_handle=output_handle,
            )
            retval["status"] = retval["status"] and ok
        else:
            if verbose:
                print(f"  Applying resize({new_shape})", file=output_handle)
            with tiledbsoma.SparseNDArray.open(
                item.uri, "w", context=context
            ) as writer:
                writer.resize(new_shape)

    elif isinstance(item, tiledbsoma.DenseNDArray):

        _print_leaf_node_banner(
            uri=item.uri,
            type_name="DenseNDArray",
            node_name=node_name,
            ms_name=ms_name,
            coll_name=coll_name,
            verbose=verbose,
            check_only=check_only,
            output_handle=output_handle,
        )
        if verbose or check_only:
            print(
                "  No action at this time, pending new-shape support for dense arrays in core 2.27",
                file=output_handle,
            )

    return retval


def _print_leaf_node_banner(
    *,
    uri: str,
    type_name: str,
    node_name: str,
    ms_name: str | None = None,
    coll_name: str | None = None,
    verbose: bool,
    check_only: bool,
    output_handle: Printable | None,
) -> None:
    if verbose or check_only:
        print("", file=output_handle)
        print(
            _get_leaf_node_description(
                type_name=type_name,
                node_name=node_name,
                ms_name=ms_name,
                coll_name=coll_name,
            ),
            file=output_handle,
        )
        print(f"  URI {uri}", file=output_handle)


def _get_leaf_node_description(
    type_name: str,
    node_name: str,
    *,
    ms_name: str | None,
    coll_name: str | None,
) -> str:
    # Return things like "ms/RNA/X/data"
    pieces = []
    if ms_name is not None:
        pieces.append("ms")
        pieces.append(ms_name)
    if coll_name is not None:
        pieces.append(coll_name)
    pieces.append(node_name)
    return f"[{type_name}] {'/'.join(pieces)} "


def _bannerize(
    node_name: str,
    value: Any,
    *,
    verbose: bool,
    check_only: bool,
    output_handle: Printable | None,
) -> None:
    if verbose or check_only:
        print(f"  {node_name:<20} {value}", file=output_handle)


def _print_dry_run_result(
    *,
    ok: bool,
    msg: str,
    output_handle: Printable | None,
) -> None:
    if ok:
        print("  OK", file=output_handle)
    else:
        print(f"  Not OK: {msg}", file=output_handle)


def _get_new_var_shape(
    *,
    nvars: dict[str, int] | None,
    ms_name: str | None,
) -> int:
    """Maps from experiment-level nobs and per-measurement nvar to the correct resizes
    for var dataframes within a SOMA experiment.
    """
    if nvars is None:
        raise tiledbsoma.SOMAError("experiment resize: internal error: nvars missing")

    if ms_name is None:
        raise tiledbsoma.SOMAError("experiment resize: internal error: ms_name missing")

    if ms_name not in nvars:
        raise tiledbsoma.SOMAError(
            f"experiment resize: missing measurement name '{ms_name}' in provided nvars"
        )
    return nvars[ms_name]


def _get_new_ndarray_shape(
    *,
    current_shape: tuple[int, ...],
    nobs: int | None,
    nvars: dict[str, int] | None,
    ms_name: str | None,
    coll_name: str | None,
) -> tuple[int, ...]:
    """Maps from experiment-level nobs and per-measurement nvar to the correct resizes
    for various kinds of n-d array within a SOMA Experiment.
    """
    if nobs is None:
        raise tiledbsoma.SOMAError("experiment resize: internal error: nobs missing")

    if nvars is None:
        raise tiledbsoma.SOMAError("experiment resize: internal error: nvars missing")

    if ms_name is None:
        raise tiledbsoma.SOMAError("experiment resize: internal error: ms_name missing")

    if ms_name not in nvars:
        raise tiledbsoma.SOMAError(
            f"experiment resize: missing measurement name '{ms_name}' in provided nvars"
        )
    nvar = nvars[ms_name]

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


def _leaf_visitor_get_shapes(
    item: Any,
    *,
    node_name: str,
    nobs: int | None,
    nvars: dict[str, int] | None,
    ms_name: str | None,
    coll_name: str | None,
    verbose: bool,
    check_only: bool,
    context: tiledbsoma.SOMATileDBContext | None,
    output_handle: Printable | None,
) -> dict[str, Any]:
    retval: dict[str, Any] = {}
    if isinstance(item, tiledbsoma.DataFrame):
        retval["uri"] = item.uri
        retval["type"] = "DataFrame"
        retval["count"] = item.count
        retval["non_empty_domain"] = item.non_empty_domain()
        retval["domain"] = item.domain
        retval["maxdomain"] = item.maxdomain
        retval["upgraded"] = item.tiledbsoma_has_upgraded_domain

    elif isinstance(item, tiledbsoma.SparseNDArray):
        retval["uri"] = item.uri
        retval["type"] = "SparseNDArray"
        retval["non_empty_domain"] = item.non_empty_domain()
        retval["shape"] = item.shape
        retval["maxshape"] = item.maxshape
        retval["upgraded"] = item.tiledbsoma_has_upgraded_shape

    elif isinstance(item, tiledbsoma.DenseNDArray):
        retval["uri"] = item.uri
        retval["type"] = "DenseNDArray"
        retval["shape"] = item.shape
        retval["maxshape"] = item.maxshape
        retval["upgraded"] = item.tiledbsoma_has_upgraded_shape

    return retval


# Expected input is like
# {
#   "obs": { "status": true },
#   "ms": {
#     "RNA": {
#       "var": { "status": true },
#       "X": {
#         "data": { "status": true }
#       }
#     }
#   }
# }
def _check_statuses(dikt: dict[str, Any]) -> bool:
    """This reduces the pass/fail statuses for all leaf nodes in a treewalk
    over the experiment down to a single pass/fail. It returns True
    only when all leaves have status=True, else False.
    """
    if "status" in dikt:
        ok = dikt["status"]
        assert isinstance(ok, bool)
        return ok
    for key, value in dikt.items():
        assert isinstance(value, dict)
        if not _check_statuses(value):
            return False
    return True
