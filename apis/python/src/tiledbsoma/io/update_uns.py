# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

from __future__ import annotations

from os.path import join
from typing import Literal

import numpy as np
import pandas as pd
from somacore.options import PlatformConfig

from tiledbsoma import Experiment, SOMATileDBContext
from tiledbsoma._collection import AnyTileDBCollection, Collection
from tiledbsoma.io._common import AdditionalMetadata, UnsMapping
from tiledbsoma.io._registration import AxisIDMapping
from tiledbsoma.io.ingest import (
    IngestionParams,
    IngestPlatformCtx,
    _ingest_uns_array,
    _maybe_set,
    _write_dataframe,
)
from tiledbsoma.logging import logger

Strict = Literal[True, "raise", "warn", "info", "debug", "dry_run"]


def update_uns_by_uri(
    uri: str,
    uns: UnsMapping,
    measurement_name: str,
    *,
    use_relative_uri: bool | None = None,
    context: SOMATileDBContext | None = None,
    additional_metadata: AdditionalMetadata = None,
    platform_config: PlatformConfig | None = None,
    default_index_name: str | None = None,
    strict: Strict = True,
) -> None:
    """Wrapper around :func:`_update_uns` that opens the experiment at the given URI.

    Some "update uns" operations, that we may want to support in the future, will require opening
    the ``Experiment`` more than once (e.g. to modify a DataFrame or array in ways that can't be
    achieved "in-place" / at one TileDB timestamp). This function API is a sketch of what that
    might look like, though it just calls :func:`update_uns` for now (overwriting DataFrames and
    arrays is currently not supported, and will either ``raise`` or "INFO" log, based on the value
    of the ``strict`` arg).
    """
    with Experiment.open(uri, "w") as exp:
        _update_uns(
            exp,
            uns,
            measurement_name,
            use_relative_uri=use_relative_uri,
            context=context,
            additional_metadata=additional_metadata,
            platform_config=platform_config,
            default_index_name=default_index_name,
            strict=strict,
        )


def _update_uns(
    exp: Experiment,
    uns: UnsMapping,
    measurement_name: str,
    *,
    use_relative_uri: bool | None = None,
    context: SOMATileDBContext | None = None,
    additional_metadata: AdditionalMetadata = None,
    platform_config: PlatformConfig | None = None,
    default_index_name: str | None = None,
    strict: Strict = True,
) -> None:
    """
    Update the given experiment/measurement's ``uns`` Collection.

    ``uns`` (short for unstructured data) is conceptually a dictionary; values can be scalars,
    DataFrames, arrays (dense or sparse), or recursively-nested dictionaries.

    This function makes a best effort at changing the existing ``uns`` Collection to match the
    provided ``uns`` dictionary. It refuses to update DataFrame and array nodes, with errors
    either raised or logged, based on the ``strict`` argument.

    Args:
        exp: The :class:`SOMAExperiment` whose ``uns`` is to be updated. Must
        be opened for write.

        measurement_name: Specifies which measurement's ``uns`` within the experiment
        is to be updated.

        uns: a Pandas dataframe with the desired contents.

        use_relative_uri: If True, store the URI relative to the experiment's URI.

        context: Optional :class:`SOMATileDBContext` containing storage parameters, etc.

        additional_metadata: Additional metadata to be added to the collection.

        platform_config: Platform-specific options used to update this array, provided
        in the form ``{"tiledb": {"create": {"dataframe_dim_zstd_level": 7}}}``

        default_index_name: Fallback name to use for columns representing ``pd.DataFrame`` indices.

        strict: How to handle conflicts with existing nodes. By default (``strict=True``), a first
        pass checks for conflicts, ``raise``ing if any are found. If none are found, a second pass
        then performs the updates. This is equivalent to running once with ``strict="dry_run"``,
        then again with ``dry_run="raise"`` (though the expectation is the latter will never
        actually ``raise``). If ``strict={debug,info,warn}`` instead, conflicting nodes are logged
        at the corresponding level, and the update is skipped.
    """
    if measurement_name not in exp.ms:
        raise ValueError(
            f"cannot find measurement name {measurement_name} within experiment at {exp.uri}"
        )

    ingest_platform_ctx = IngestPlatformCtx(
        context=context,
        ingestion_params=IngestionParams("write", label_mapping=None),
        additional_metadata=additional_metadata,
        platform_config=platform_config,
    )
    if strict is True:
        # Surface any errors (e.g. unsupported SOMA object overwrites) without writing anything
        _update_uns_dict(
            exp.ms[measurement_name]["uns"],  # type: ignore[arg-type]
            uns,
            use_relative_uri=use_relative_uri,
            ingest_platform_ctx=ingest_platform_ctx,
            default_index_name=default_index_name,
            strict="dry_run",
        )
        strict = "raise"
    _update_uns_dict(
        exp.ms[measurement_name]["uns"],  # type: ignore[arg-type]
        uns,
        use_relative_uri=use_relative_uri,
        ingest_platform_ctx=ingest_platform_ctx,
        default_index_name=default_index_name,
        strict=strict,
    )


def _update_uns_dict(
    coll: AnyTileDBCollection,
    uns: UnsMapping,
    *,
    ingest_platform_ctx: IngestPlatformCtx,
    use_relative_uri: bool | None,
    default_index_name: str | None,
    strict: Strict,
) -> None:
    for k, v in uns.items():
        # Any existing scalar will be found in metadata, not as a child of the Collection.
        cur = None
        if k in coll:
            cur = coll[k]
        if k in coll.metadata:
            if cur is not None:
                logger.warn(f"{coll.uri}[{k}] exists as both metadata and child")
            else:
                cur = coll.metadata[k]
        exists = cur is not None

        def can_write() -> bool:
            if exists:
                msg = f"{coll.uri}[{k}]: already exists (type {type(cur).__name__}), refusing to overwrite with {v}"
                if strict in ["dry_run", "raise"]:
                    raise ValueError(msg)
                else:
                    msg = f"Skipping {msg}"
                    if strict == "warn":
                        logger.warn(msg)
                    elif strict == "info":
                        logger.info(msg)
                    elif strict == "debug":
                        logger.debug(msg)
                    return False
            else:
                return strict != "dry_run"

        if isinstance(v, (str, int, float, np.generic)):
            if k in coll:
                raise ValueError(
                    f"can't overwrite {type(cur).__name__} at {coll.uri}/{k} with scalar {v}"
                )
            if isinstance(v, np.generic):
                # Unwrap numpy scalar
                v = v.item()
            if strict != "dry_run":
                coll.metadata[k] = v
        elif isinstance(v, pd.DataFrame):
            if can_write():
                with _write_dataframe(
                    df_uri=join(coll.uri, k),
                    df=v.copy(),  # `_write_dataframe` modifies the `pd.DataFrame` it's passed
                    id_column_name=default_index_name,
                    ingestion_params=ingest_platform_ctx["ingestion_params"],
                    additional_metadata=ingest_platform_ctx["additional_metadata"],
                    platform_config=ingest_platform_ctx["platform_config"],
                    context=ingest_platform_ctx["context"],
                    axis_mapping=AxisIDMapping.identity(v.shape[0]),
                ) as df:
                    _maybe_set(coll, k, df, use_relative_uri=use_relative_uri)
        elif isinstance(v, dict):
            if exists:
                if not isinstance(cur, Collection):
                    raise ValueError(
                        f"{coll.uri}/{k}: expected Collection, found {type(cur).__name__}"
                    )
            _update_uns_dict(
                coll[k],
                v,
                use_relative_uri=use_relative_uri,
                ingest_platform_ctx=ingest_platform_ctx,
                default_index_name=default_index_name,
                strict=strict,
            )
        elif isinstance(v, np.ndarray):
            if can_write():
                _ingest_uns_array(
                    coll,
                    k,
                    v,
                    use_relative_uri=use_relative_uri,
                    ingest_platform_ctx=ingest_platform_ctx,
                )
        else:
            raise ValueError(f"unsupported uns type {type(v)}")
