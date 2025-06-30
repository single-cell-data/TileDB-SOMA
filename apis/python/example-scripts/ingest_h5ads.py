#!/usr/bin/env python

"""
Demo concurrent ingestion of H5AD using multi-processing (concurrent.futures.ProcessPoolExecutor).

See also the Academy bulk-ingestion tutorial here:
https://documentation.cloud.tiledb.com/academy/structure/life-sciences/single-cell/tutorials/

Demo assumptions:
* All H5AD can be opened with `anndata.read_h5ad`, i.e., they are available in a file system
* The newly created SOMA Experiment can be written to any TileDB-supported storage system (e.g., S3, etc).

There are two workflows demonstrated -- the `--normalize` command-line argument controls which
is used:

    1. --no-normalize: Load data directly from H5ADs. This workflow is slightly faster
       when the H5ADs can be loaded "as is" in the registration step.

    2. --normalize: Load from AnnData -- in this case, the workflow includes a step where
       user-specified code loads and modifies the AnnData. While this approach is slower,
       it is appropriate when the workflow requires on-the-fly AnnData schema normalization
       or other modifications.

The first workflow runs by default.

This is just a demo -- it omits features such as loading additional X layers, etc. It should
port without substantial changes to other concurrent computing frameworks.

Example usage:

    $ python ingest_h5ads.py --experiment-uri [path] --h5ad-directory [path] --num-workers 8

Notes on usage:

1. Memory use will increase with more workers. Total memory uses is determined by both `num-workers`
   and the memory required to load and process each AnnDAta (each worker loads a separate AnnData).
   Each AnnData is loaded in "backed" mode, which includes the
   `obs` and `var` dataframes, the `obsm`/`obsp`/`varm`/`varp` matrices, as well as the `uns` dictionary.

   In other words, the total memory budget is the AnnData memory usage multiplied by the number of
   workers.

   Every data collection is different, and some experimentation will be required to determine the
   memory budget for each worker, and from that the optimal number of workers.  As a place to start,
   it is common that large H5ADs, _without_ obsm/varm/obsp/varp matrices, will require roughly 16-20GiB
   to process, and smaller files less. Put another way, on a 64GiB host, processing large H5ADs,
   you might want to start with 3 or 4 workers.  If your datasets are smaller, this could be increased.

2. While `anndata.read_h5ad()` is unable to read data directly from S3, the `tiledbsoma.io.from_h5ad()` and
   `tiledbsoma.io.register_h5ad()` functions will accept any URL supported by TileDB VFS (including S3).

"""

from __future__ import annotations

import argparse
import concurrent.futures
import gc
import logging
import logging.config
import multiprocessing
import os
import sys
from concurrent.futures import Executor, Future, ProcessPoolExecutor
from itertools import repeat
from typing import Any, Callable, Iterable, Iterator, TypeVar

import anndata as ad

import tiledbsoma
import tiledbsoma.io

try:
    from tiledb.cloud import Config as CloudConfig

    TileDBCloudConfig = CloudConfig().dict()

except (ImportError, ModuleNotFoundError):
    TileDBCloudConfig = {}

logger = logging.getLogger("ingest_h5ads")


def main():
    # "fork" is inherently unsafe on Posix/Linux systems when both threads and sub-processes are used.
    # Use "spawn", which is safe, and will be the default on all systems in the near future.
    multiprocessing.set_start_method("spawn", force=True)

    logger_setup()

    args = parse_args()
    logger.info(repr(args))

    h5ad_paths = get_h5ad_paths(args.h5ad_directory)
    logger.info(f"Found {len(h5ad_paths)} H5ADs to ingest")

    experiment_uri = args.experiment_uri
    # work-around sc-64882 - experiment URI must not have a trailing slash.
    if experiment_uri[-1] == "/":
        experiment_uri = experiment_uri[0:-1]

    context = tiledbsoma.SOMATileDBContext(tiledb_config=TileDBCloudConfig)

    ##
    # Step 1 - create Experiment (empty)
    ##
    if not args.append:
        tiledbsoma.io.from_anndata(
            experiment_uri,
            (
                ad.read_h5ad(h5ad_paths[0], backed="r")
                if not args.normalize
                else normalize_anndata(ad.read_h5ad(h5ad_paths[0], backed="r"))
            ),
            measurement_name=args.measurement_name,
            obs_id_name=args.obs_id_name,
            var_id_name=args.var_id_name,
            X_layer_name=args.X_layer_name,
            ingest_mode="schema_only",
            context=context,
        )
        logger.info("Experiment created")

    ##
    # Step 2 - build ingestion plan
    ##
    logger.info("Registering all H5ADs")
    if args.normalize:
        registration_mapping = tiledbsoma.io.register_anndatas(
            experiment_uri,
            (normalize_anndata(ad.read_h5ad(p, backed="r")) for p in h5ad_paths),
            measurement_name=args.measurement_name,
            obs_field_name=args.obs_id_name,
            var_field_name=args.var_id_name,
            context=context,
        )

    else:
        registration_mapping = tiledbsoma.io.register_h5ads(
            experiment_uri,
            h5ad_paths,
            measurement_name=args.measurement_name,
            obs_field_name=args.obs_id_name,
            var_field_name=args.var_id_name,
            context=context,
            use_multiprocessing=True,  # performance improvement when reading H5AD files
        )
    logger.info(
        f"Ingest plan created, n_obs={registration_mapping.get_obs_shape()}, n_vars={registration_mapping.get_var_shapes()}",
    )

    ##
    # Step 3 - evolve the experiment in preparation for parallel ingest
    ##
    registration_mapping.prepare_experiment(experiment_uri, context=context)
    logger.info("Experiment prepared for ingestion.")

    ##
    # Step 4 - parallel ingest each H5AD
    ##
    with ProcessPoolExecutor(max_workers=args.num_workers, initializer=logger_setup, max_tasks_per_child=1) as ppe:
        # Use a lazy executor so that the params which consume a lot of memory (e.g.,registration_map)
        # are materialized as late as is possible.
        #
        # In addition, call `subset_for_h5ad` to reduce the size of the object transmitted
        # to the distributed worker.
        for result in LazyExecutor(ppe, args.num_workers + max(1, args.num_workers // 2)).map(
            worker_load_dataset,
            h5ad_paths,
            repeat(context.tiledb_config),
            repeat(experiment_uri),
            repeat(args.measurement_name),
            repeat(args.obs_id_name),
            repeat(args.var_id_name),
            repeat(args.X_layer_name),
            (registration_mapping.subset_for_h5ad(fn) for fn in h5ad_paths),  # lazy
            repeat(args.normalize),
            ((idx + 1, len(h5ad_paths)) for idx in range(len(h5ad_paths))),
        ):
            gc.collect()
            logger.debug(f"Completed {result}")

    logger.info("Ingestion complete")
    return 0


def normalize_anndata(adata: ad.AnnData) -> ad.AnnData:
    """Hook to perform any required AnnData on-the-fly "normalization" steps, at dataloading
    time, such as removing obsp/varp, cleaning up uns, etc.

    This code is centralized to ensure we perform the same normalization steps on all ingestion steps,
    and will be called during Experiment creation (schema determination), AnnData registration, and
    per-H5AD ingestion.

    This is only used if `--normalize` option is set.
    """

    # Customize as needed. By default, removing all AnnData slots that do not naturally
    # support simple concatenation.
    del adata.obsm
    del adata.varm
    del adata.obsp
    del adata.varp
    del adata.uns
    del adata.layers

    return adata


def worker_load_dataset(
    filename: str,
    tiledb_config: dict[str, Any],
    experiment_uri,
    measurement_name: str,
    obs_id_name: str,
    var_id_name: str,
    X_layer_name: str,
    registration_mapping,
    normalize: bool,
    i_of_j: tuple[int, int],
) -> tuple[tuple[int, int], str]:
    """Worker callback to ingest a single H5AD.

    Returns (i_of_j, filename).
    """
    i, j = i_of_j
    logger.info(f"Worker start [{i} of {j}] {filename}")
    context = tiledbsoma.SOMATileDBContext(tiledb_config=tiledb_config)

    if normalize:
        # load AnnData and call the user-defined normalized hook
        adata = normalize_anndata(ad.read_h5ad(filename))
        tiledbsoma.io.from_anndata(
            experiment_uri,
            adata,
            measurement_name=measurement_name,
            obs_id_name=obs_id_name,
            var_id_name=var_id_name,
            X_layer_name=X_layer_name,
            uns_keys=(),
            registration_mapping=registration_mapping,
            context=context,
        )
        del adata

    else:
        # just ingest the H5AD directly
        tiledbsoma.io.from_h5ad(
            experiment_uri,
            filename,
            measurement_name=measurement_name,
            obs_id_name=obs_id_name,
            var_id_name=var_id_name,
            X_layer_name=X_layer_name,
            uns_keys=(),
            registration_mapping=registration_mapping,
            context=context,
        )

    del context
    gc.collect()
    logger.info(f"Worker done [{i} of {j}] {filename}")
    return i_of_j, filename


_T = TypeVar("_T")


class LazyExecutor:
    """This class provides the functionality of `concurrent.futures.Executor.map`, but
    extracts arguments lazily rather than eagerly. This is important for pipelines that
    need to minimize memory use while doing very large concurrent mappings.

    If you currently do:

        >>> results = executor.map(fn, args)

    replace with:

        >>> results = LazyExecutorMap(executor, num_pending).map(fn, args)

    where `num_pending` is the max number of submitted tasks at any given time. Typically
    `num_pending` should have the same value as the Executor `max_workers` parameter, or a
    larger value if queued work is desired.

    Works with both ThreadPoolExecutor and ProcessPoolExecutor
    """

    def __init__(self, executor: Executor, num_pending: int) -> None:
        self.ex = executor
        self.num_pending = num_pending

    def map(self, fn: Callable[..., _T], *iterables: Iterable[Any]) -> Iterator[_T]:
        itr = iter(zip(*iterables))
        pending: list[Future[_T]] = []
        in_progress: set[Future[_T]] = set()
        try:
            while True:
                for _ in range(self.num_pending - len(in_progress)):
                    el = next(itr)
                    fut = self.ex.submit(fn, *el)
                    pending.append(fut)
                    in_progress.add(fut)
                    logger.debug(f"Queued job, num pending: {len(pending)}")

                _, in_progress = concurrent.futures.wait(in_progress, return_when=concurrent.futures.FIRST_COMPLETED)

                while pending and pending[0].done():
                    yield pending.pop(0).result()

        except StopIteration:
            in_progress.clear()
            while pending:
                yield pending.pop(0).result()


def cpu_count() -> int:
    return os.cpu_count() or 1


def default_cpu_count() -> int:
    return max(cpu_count() // 4, 1)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment-uri", type=str, help="Experiment URI", required=True)
    parser.add_argument("--h5ad-directory", type=str, help="Folder of H5AD files", required=True)
    parser.add_argument(
        "--append",
        default=False,
        action=argparse.BooleanOptionalAction,
        help="Append to existing Experiment (default: no-append)",
    )
    parser.add_argument(
        "--normalize",
        default=False,
        action=argparse.BooleanOptionalAction,
        help="Demonstrate pre-ingest normalization of AnnData (default: off)",
    )
    parser.add_argument(
        "--num-workers",
        type=int,
        default=0,
        help=f"Number of ingestion workers - if zero, will use all available cores (default: {default_cpu_count()})",
    )
    parser.add_argument(
        "--measurement-name",
        type=str,
        default="RNA",
        help="Measurement name used during ingest",
    )
    parser.add_argument(
        "--obs-id-name",
        type=str,
        default="obs_id",
        help="obs column name containing unique identifiers",
    )
    parser.add_argument(
        "--var-id-name",
        type=str,
        default="var_id",
        help="var column name containing unique identifiers",
    )
    parser.add_argument("--X-layer-name", type=str, default="data", help="Default X layer name")

    args = parser.parse_args()

    if args.num_workers <= 0:
        args.num_workers = cpu_count()
    if args.num_workers > cpu_count():
        args.num_workers = cpu_count()

    return args


def logger_setup() -> None:
    """set up the Python logger"""
    logging.config.dictConfig(
        {
            "version": 1,
            "formatters": {
                "default": {
                    "class": "logging.Formatter",
                    "format": "[%(asctime)s %(process)d] %(message)s",
                },
            },
            "handlers": {"console": {"class": "logging.StreamHandler", "formatter": "default"}},
            "loggers": {
                "ingest_h5ads": {"level": "INFO", "handlers": ["console"]},
                "tiledbsoma": {"level": "WARNING", "handlers": ["console"]},
            },
        },
    )


def get_h5ad_paths(h5ad_directory_path: str) -> list[str]:
    """Return list of all H5AD paths in the user-provided folder."""
    return [
        os.path.join(h5ad_directory_path, fn)
        for fn in os.listdir(h5ad_directory_path)
        if fn.endswith(".h5ad") and os.path.isfile(os.path.join(h5ad_directory_path, fn))
    ]


if __name__ == "__main__":
    sys.exit(main())
