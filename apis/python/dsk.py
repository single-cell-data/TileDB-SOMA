#!/usr/bin/env python

import json
from contextlib import nullcontext
from functools import wraps
from os import cpu_count, makedirs
from os.path import dirname, join, splitext
from sys import stderr, stdout

import dask
import pyarrow as pa
import scanpy as sc
from anndata import AnnData
from click import argument, group, option
from dask.distributed import Client, LocalCluster
from humanfriendly import parse_size
from prometheus_client.decorator import getfullargspec
from pyarrow import feather
from scipy.sparse import csr_matrix, vstack
from somacore import AxisQuery
from utz import Time, err, iec
from utz.cli import number
from utz.mem import Tracker

from tiledbsoma import Experiment, SparseNDArray
from tiledbsoma._dask import DaskConfig, make_context
from tiledbsoma._fastercsx import CompressedMatrix
from tiledbsoma._indexer import IntIndexer


def sz(size: int) -> str:
    return f"{size} ({iec(size)})"


@group
def cli():
    pass


CENSUS_S3 = "s3://cellxgene-census-public-us-west-2/cell-census"
CENSUS = f"{CENSUS_S3}/2024-07-01/soma/census_data/homo_sapiens"


def _census_version_opt(*opt_args, **opt_kwargs):
    def opt(fn):
        @option(*opt_args, "--census-version", "census_version", **opt_kwargs)
        @wraps(fn)
        def _fn(*args, census_version, **kwargs):
            exp_uri = f"{CENSUS_S3}/{census_version}/soma/census_data/homo_sapiens"
            return fn(*args, exp_uri=exp_uri, **kwargs)

        return _fn

    return opt


# fmt: off
mem_budget_opt = option("-b", "--mem-total-budget", callback=lambda ctx, param, val: parse_size(val) if val else None)
dask_chunk_size_opt = number("-c", "--dask-chunk-size", default="100k")
census_version_opt = _census_version_opt("-C", default="2024-07-01")
no_dask_mem_mgmt_opt = option("-D", "--no-dask-mem-mgmt", is_flag=True)
no_keep_memray_bin_opt = option("-K", "--no-keep-memray-bin", is_flag=True)
x_layer_opt = option("-l", "--x-layer", default="raw")
memray_bin_opt = option("-m", "--memray-bin-path")
memray_level_opt = option("-M", "--memray-level", count=True, help="0x: no memray, 1x: memray, 2x: memray with `trace_python_allocators` (slower)")
method_opt = option("-M", "--method", count=True, help='0x: "naive" reindexing+CSR-construction, 1x: .blockwise().scipy, 2x: IntIndexer/fastercsx')
out_dir_opt = option("-o", "--out-dir")
dask_dashboard_port_opt = option("-p", "--dask-dashboard-port", type=int, default=8787)
tissue_opt = option("-t", "--tissue", default="nose")
tdb_concurrency_opt = option("-T", "--tdb-concurrency", type=int, default=1)
verbosity_opt = option("-v", "--verbosity", default=2)
var_filter_opt = option("-V", "--var-filter")
dask_workers_opt = option("-w", "--n-dask-workers", type=int)
dask_threads_per_worker_opt = option("-W", "--n-dask-threads-per-worker", type=int)

hvg_inplace_opt = option("--hvg-inplace/--no-hvg-inplace", default=True)
hvg_subset_opt = option("--hvg-subset/--no-hvg-subset", default=False)
# fmt: on


DEFAULT_CONFIG = {
    "vfs.s3.no_sign_request": "true",
    "vfs.s3.region": "us-west-2",
}


@cli.command
@tissue_opt
@tdb_concurrency_opt
@argument("joinids_dir")
def joinids(
    tissue: str,
    tdb_concurrency: int | None,
    joinids_dir: str,
):
    """Output feather files containing the obs and var joinids responsive to a CELLxGENE Census query."""
    time = Time(log=True)
    makedirs(joinids_dir, exist_ok=True)
    obs_path = join(joinids_dir, "obs.feather")
    var_path = join(joinids_dir, "var.feather")
    context = make_context(
        tdb_concurrency=tdb_concurrency,
        tiledb_config=DEFAULT_CONFIG,
    )
    time("open")
    with Experiment.open(CENSUS, "r", context=context) as exp:
        time("query")
        query = exp.axis_query(
            measurement_name="RNA",
            obs_query=AxisQuery(
                value_filter=f'is_primary_data == True and tissue_general == "{tissue}"'
            ),
        )
        time("obs_joinids")
        obs_joinids = query.obs_joinids()
        time("obs_tbl")
        obs_tbl = pa.Table.from_arrays([obs_joinids], names=["obs_joinids"])
        time("obs_save")
        feather.write_feather(obs_tbl, obs_path)
        time()
        err(f"Wrote {len(obs_joinids):,} obs joinids to {obs_path}")
        time("var_joinids")
        var_joinids = query.var_joinids()
        time("var_tbl")
        var_tbl = pa.Table.from_arrays([var_joinids], names=["var_joinids"])
        time("var_save")
        feather.write_feather(var_tbl, var_path)
        time()
        err(f"Wrote {len(var_joinids):,} var joinids to {var_path}")

        json.dump(
            time.fmt(),
            stderr,
            indent=2,
        )
        err()


def memray_cmd(out_path_kwarg: str):
    """Decorator for dsk.py subcommands that instrument Memray profiling."""

    def rv(fn):
        @cli.command
        @dask_chunk_size_opt
        @no_keep_memray_bin_opt
        @memray_bin_opt
        @memray_level_opt
        @out_dir_opt
        @tdb_concurrency_opt
        @verbosity_opt
        @wraps(fn)
        def _fn(
            *args,
            dask_chunk_size: int,
            no_keep_memray_bin: bool,
            memray_bin_path: str | None,
            memray_level: int,
            out_dir: str | None,
            tdb_concurrency: int,
            verbosity: int,
            **kwargs,
        ):
            if memray_bin_path is None:
                if dask_chunk_size % 1000 == 0:
                    name = f"{dask_chunk_size // 1000}k"
                else:
                    name = f"{dask_chunk_size}"

                if tdb_concurrency == 0:
                    tdb_concurrency = cpu_count()

                if "n_dask_workers" in kwargs and "n_dask_threads_per_worker" in kwargs:
                    n_dask_workers = kwargs["n_dask_workers"]
                    n_dask_threads_per_worker = kwargs["n_dask_threads_per_worker"]
                    name += f"_{n_dask_workers}x{n_dask_threads_per_worker}x{tdb_concurrency}"
                elif tdb_concurrency != 1:
                    name += f"_{tdb_concurrency}"
                memray_bin_path = join(out_dir or kwargs[out_path_kwarg], f"{name}.bin")

            out_stem = splitext(memray_bin_path)[0]

            time = Time()

            if not tdb_concurrency:
                tdb_concurrency = cpu_count()

            spec = getfullargspec(fn)
            wrapped = getattr(fn, "__wrapped__", None)
            wrapped_spec = getfullargspec(wrapped) if wrapped else None
            arg_ks = spec.args + (wrapped_spec.args if wrapped_spec else [])
            opt_kwargs = dict(
                dask_chunk_size=dask_chunk_size,
                out_dir=out_dir,
                tdb_concurrency=tdb_concurrency,
                out_stem=out_stem,
            )

            if not memray_level:
                mem = None
                mem_ctx = nullcontext()
            else:
                opt_kwargs["memray_bin_path"] = memray_bin_path
                mem = mem_ctx = Tracker(
                    memray_bin_path,
                    native_traces=True,
                    follow_fork=True,
                    trace_python_allocators=memray_level > 1,
                    keep=not no_keep_memray_bin,
                    log=verbosity,
                )

            opt_kwargs = {k: v for k, v in opt_kwargs.items() if k in arg_ks}

            time("start")
            return fn(
                *args,
                **opt_kwargs,
                time=time,
                mem=mem,
                mem_ctx=mem_ctx,
                **kwargs,
            )

    return rv


def query_to_anndata_dask(
    uri: str,
    tiledb_config: dict[str, str],
    mem_total_budget: str | None,
    dask_chunk_size: int,
    x_layer: str,
    dask_dashboard_port: int,
    time: Time,
    tissue: str,
    tdb_concurrency: int,
    var_filter: str | None,
    no_dask_mem_mgmt: bool,
    n_dask_workers: int | None,
    n_dask_threads_per_worker: int | None,
    normalize: bool = True,
    log1p: bool = True,
) -> AnnData:
    if mem_total_budget:
        tiledb_config.update(**{"sm.mem.total_budget": mem_total_budget})
    context = make_context(
        tiledb_config=tiledb_config,
        tdb_concurrency=tdb_concurrency,
    )
    exp = Experiment.open(
        uri,
        context=context,
    )
    obs_filter = "is_primary_data == True"
    if tissue:
        obs_filter = f'tissue_general == "{tissue}" and {obs_filter}'

    time("dask")
    if dask_chunk_size:
        if no_dask_mem_mgmt:
            dask.config.set(
                {
                    # "distributed.worker.memory.target": 0.7,
                    "distributed.worker.memory.spill": False,
                    "distributed.worker.memory.pause": False,
                    "distributed.worker.memory.terminate": False,
                }
            )
        if n_dask_workers is None:
            n_dask_workers = cpu_count() // tdb_concurrency // n_dask_threads_per_worker
        cluster = LocalCluster(
            n_workers=n_dask_workers,
            threads_per_worker=n_dask_threads_per_worker,
            dashboard_address=f":{dask_dashboard_port}",
        )
        err(f"{n_dask_workers=}, {n_dask_threads_per_worker=}, {tdb_concurrency=}")
        client = Client(cluster)
        err(f"Dask client: {client}, dashboard: {client.dashboard_link}")
    else:
        dask_chunk_size = None

    err("Running query")
    time("query")
    query = exp.axis_query(
        measurement_name="RNA",
        obs_query=AxisQuery(value_filter=obs_filter) if obs_filter else None,
        var_query=AxisQuery(value_filter=var_filter) if var_filter else None,
    )
    time("to_anndata")
    add = query.to_anndata(
        x_layer,
        dask=DaskConfig(
            chunk_size=dask_chunk_size,
            tdb_concurrency=tdb_concurrency,
            tdb_configs=tiledb_config,
        ),
    )
    if normalize:
        err("normalize...")
        time("normalize")
        sc.pp.normalize_total(add)
    if log1p:
        err("log1p...")
        time("log1p")
        sc.pp.log1p(add)
    time()
    return add


def process_mem(
    key: str,
    tissue: str,
    mem: Tracker | None,
    time: Time,
    out_stem: str | None,
    n_obs: int,
    mem_total_budget: str | None,
    no_dask_mem_mgmt: bool,
    dask_chunk_size: int,
    n_dask_workers: int,
    n_dask_threads_per_worker: int,
    tdb_concurrency: int,
):
    mem_kwargs = {}
    if mem:
        peak_mem = mem.peak_mem
        mem_kwargs["peak_mem"] = peak_mem
        err(f"Peak memory use: {sz(peak_mem)}")

        # Copy some memray stats, include IEC string reprs (e.g. "2.1 GiB")
        stats = mem.stats
        top_allocs = [
            {**alloc, "iec": iec(alloc["size"])}
            for alloc in stats["top_allocations_by_size"]
        ]
        mem_kwargs.update(
            {
                "top_allocs": top_allocs,
                "tracked_allocs": stats["total_num_allocations"],
                "tracked_bytes": stats["total_bytes_allocated"],
                "total_allocs": stats["metadata"]["total_allocations"],
                "total_frames": stats["metadata"]["total_frames"],
            }
        )

    total_measured_time = sum(time.times.values())
    obs_per_s = n_obs / total_measured_time
    key_obs_per_s = n_obs / time[key]
    key_obs_per_s_key = f"{key}_obs_per_s"
    err(
        f"Total measured time: {total_measured_time:.4g} ({round(obs_per_s):,} cells/s, {key.upper()}: {round(key_obs_per_s):,} cells/s)"
    )
    out_json_path = f"{out_stem}.json"
    stats = dict(
        **mem_kwargs,
        times=time.times,
        total_time=total_measured_time,
        tissue=tissue,
        cmd=key.upper(),
        n_obs=n_obs,
        obs_per_s=obs_per_s,
        **({"no_dask_mem_mgmt": no_dask_mem_mgmt} if no_dask_mem_mgmt else {}),
        dask_chunk_size=dask_chunk_size,
        n_dask_workers=n_dask_workers,
        n_dask_threads_per_worker=n_dask_threads_per_worker,
        tdb_concurrency=tdb_concurrency,
        **{key_obs_per_s_key: key_obs_per_s},
        **({"mem_total_budget": mem_total_budget} if mem_total_budget else {}),
    )
    json.dump(stats, stdout, indent=2)
    print()
    makedirs(dirname(out_json_path), exist_ok=True)
    with open(out_json_path, "w") as f:
        json.dump(stats, f, indent=2)
        print(file=f)


@memray_cmd(out_path_kwarg="tissue")
@mem_budget_opt
@census_version_opt
@hvg_inplace_opt
@hvg_subset_opt
@x_layer_opt
@dask_dashboard_port_opt
@tissue_opt
@var_filter_opt
@no_dask_mem_mgmt_opt
@dask_workers_opt
@dask_threads_per_worker_opt
def hvg(
    mem_total_budget: str | None,
    dask_chunk_size: int,
    hvg_inplace: bool,
    hvg_subset: bool,
    x_layer: str,
    out_stem: str,
    dask_dashboard_port: int,
    exp_uri: str,
    time: Time,
    mem: Tracker | None,
    mem_ctx: Tracker | nullcontext,
    tissue: str,
    tdb_concurrency: int,
    var_filter: str | None,
    no_dask_mem_mgmt: bool,
    n_dask_workers: int | None,
    n_dask_threads_per_worker: int | None,
):
    """Run Scanpy HVG on an AnnData whose ``X`` array is a Dask-backed TileDB-SOMA CELLxGENE Census query."""
    with mem_ctx:
        add = query_to_anndata_dask(
            uri=exp_uri,
            tiledb_config=DEFAULT_CONFIG,
            mem_total_budget=mem_total_budget,
            dask_chunk_size=dask_chunk_size,
            x_layer=x_layer,
            dask_dashboard_port=dask_dashboard_port,
            time=time,
            tissue=tissue,
            tdb_concurrency=tdb_concurrency,
            var_filter=var_filter,
            no_dask_mem_mgmt=no_dask_mem_mgmt,
            n_dask_workers=n_dask_workers,
            n_dask_threads_per_worker=n_dask_threads_per_worker,
        )
        n_obs = add.shape[0]
        err("Running Scanpy HVG")
        time("hvg")
        sc.pp.highly_variable_genes(add, inplace=hvg_inplace, subset=hvg_subset)
        time("stats")
    time()

    process_mem(
        key="hvg",
        tissue=tissue,
        mem=mem,
        time=time,
        out_stem=out_stem,
        n_obs=n_obs,
        mem_total_budget=mem_total_budget,
        no_dask_mem_mgmt=no_dask_mem_mgmt,
        dask_chunk_size=dask_chunk_size,
        n_dask_workers=n_dask_workers,
        n_dask_threads_per_worker=n_dask_threads_per_worker,
        tdb_concurrency=tdb_concurrency,
    )


@memray_cmd(out_path_kwarg="tissue")
@mem_budget_opt
@census_version_opt
@hvg_inplace_opt
@hvg_subset_opt
@x_layer_opt
@dask_dashboard_port_opt
@tissue_opt
@var_filter_opt
@no_dask_mem_mgmt_opt
@dask_workers_opt
@dask_threads_per_worker_opt
def pca(
    mem_total_budget: str | None,
    dask_chunk_size: int,
    hvg_inplace: bool,
    hvg_subset: bool,
    x_layer: str,
    out_stem: str,
    dask_dashboard_port: int,
    exp_uri: str,
    time: Time,
    mem: Tracker | None,
    mem_ctx: Tracker | nullcontext,
    tissue: str,
    tdb_concurrency: int,
    var_filter: str | None,
    no_dask_mem_mgmt: bool,
    n_dask_workers: int | None,
    n_dask_threads_per_worker: int | None,
):
    """Run Scanpy PCA on an AnnData whose ``X`` array is a Dask-backed TileDB-SOMA CELLxGENE Census query."""
    tiledb_config = {
        "vfs.s3.no_sign_request": "true",
        "vfs.s3.region": "us-west-2",
    }
    with mem_ctx:
        add = query_to_anndata_dask(
            uri=exp_uri,
            tiledb_config=tiledb_config,
            mem_total_budget=mem_total_budget,
            dask_chunk_size=dask_chunk_size,
            x_layer=x_layer,
            dask_dashboard_port=dask_dashboard_port,
            time=time,
            tissue=tissue,
            tdb_concurrency=tdb_concurrency,
            var_filter=var_filter,
            no_dask_mem_mgmt=no_dask_mem_mgmt,
            n_dask_workers=n_dask_workers,
            n_dask_threads_per_worker=n_dask_threads_per_worker,
        )
        n_obs = add.shape[0]
        err("Running Scanpy HVG")
        time("hvg")
        sc.pp.highly_variable_genes(add, inplace=hvg_inplace, subset=hvg_subset)
        err("Running Scanpy PCA")
        time("pca")
        sc.pp.pca(add)
        time("stats")
    time()

    process_mem(
        key="pca",
        tissue=tissue,
        mem=mem,
        time=time,
        out_stem=out_stem,
        n_obs=n_obs,
        mem_total_budget=mem_total_budget,
        no_dask_mem_mgmt=no_dask_mem_mgmt,
        dask_chunk_size=dask_chunk_size,
        n_dask_workers=n_dask_workers,
        n_dask_threads_per_worker=n_dask_threads_per_worker,
        tdb_concurrency=tdb_concurrency,
    )


@memray_cmd("joinids_dir")
@mem_budget_opt
@census_version_opt
@method_opt
@argument("joinids_dir")
def csr(
    mem_total_budget: str | None,
    dask_chunk_size: int | None,
    exp_uri: str,
    memray_bin_path: str | None,
    method: int,
    time: Time,
    mem: Tracker,
    tdb_concurrency: int,
    joinids_dir: str,
):
    """Profile a simulated Dask chunk read of a slice of a CELLxGENE query."""
    with mem:
        obs_joinids_path = join(joinids_dir, "obs.feather")
        var_joinids_path = join(joinids_dir, "var.feather")
        obs_tbl = feather.read_table(obs_joinids_path)
        all_obs_joinids = obs_tbl["obs_joinids"]
        obs_joinids = all_obs_joinids[:dask_chunk_size].to_numpy().tolist()
        var_tbl = feather.read_table(var_joinids_path)
        var_joinids = var_tbl["var_joinids"].to_numpy().tolist()
        shape = (len(obs_joinids), len(var_joinids))
        soma_ctx = make_context(
            tdb_concurrency=tdb_concurrency,
            tiledb_config={
                **DEFAULT_CONFIG,
                "sm.mem.total_budget": mem_total_budget,
            },
        )
        uri = f"{exp_uri}/ms/RNA/X/raw"
        time("open")
        with SparseNDArray.open(uri, context=soma_ctx) as arr:
            if method == 0:
                time("read")
                tbl = arr.read((obs_joinids, var_joinids)).tables().concat()
                time("cols")
                soma_dim_0, soma_dim_1, data = [col.to_numpy() for col in tbl.columns]
                time("reindex-dicts")
                obs_joinid_idx_map = {
                    obs_joinid: idx
                    for idx, obs_joinid in enumerate(sorted(obs_joinids))
                }
                var_joinid_idx_map = {
                    var_joinid: idx
                    for idx, var_joinid in enumerate(sorted(var_joinids))
                }
                time("maps")
                obs = [obs_joinid_idx_map[int(obs_joinid)] for obs_joinid in soma_dim_0]
                var = [var_joinid_idx_map[int(var_joinid)] for var_joinid in soma_dim_1]
                time("csr")
                csr = csr_matrix((data, (obs, var)), shape=shape)
                time()
            elif method == 1:
                time("read")
                scipy = arr.read((obs_joinids, var_joinids)).blockwise(0).scipy()
                time("csrs")
                csrs, idxs = zip(*list(iter(scipy)))
                time("close")
                time("csr")
                csr = vstack(csrs)
                time()
                if len(csrs) > 1:
                    for i, c in enumerate(csrs):
                        err(f"CSR block {i}: {repr(c)}")
            elif method == 2:
                time("indexers")
                obs_indexer = IntIndexer(
                    obs_joinids, context=soma_ctx
                )  # TODO: i32-based reindexing?
                var_indexer = IntIndexer(var_joinids, context=soma_ctx)
                new_tbls = []
                time("read")
                tbls = arr.read((obs_joinids, var_joinids)).tables()
                idx = 0
                while True:
                    time(f"read-tbl-{idx}")
                    try:
                        tbl = next(tbls)
                    except StopIteration:
                        break
                    time(f"indexing-{idx}")
                    new_dim0 = obs_indexer.get_indexer(tbl["soma_dim_0"]).astype(
                        "int32"
                    )
                    new_dim1 = var_indexer.get_indexer(tbl["soma_dim_1"]).astype(
                        "int32"
                    )
                    time(f"new-tbl-{idx}")
                    new_tbl = pa.Table.from_pydict(
                        {
                            "soma_dim_0": new_dim0,
                            "soma_dim_1": new_dim1,
                            "soma_data": tbl["soma_data"],
                        }
                    )
                    new_tbls.append(new_tbl)
                    idx += 1
                time("CompressedMatrix")
                cm = CompressedMatrix.from_soma(
                    new_tbls,
                    shape=shape,
                    format="csr",
                    make_sorted=True,
                    context=soma_ctx,
                )
                time("csr")
                csr = cm.to_scipy()
                time()
            else:
                raise ValueError(f"Unrecognized -M/--method count: {method}")

        nnz = csr.nnz
        err(f"CSR: {csr.shape}, {csr.nnz}")

    peak = mem.peak_mem
    err(f"Peak memory use: {sz(peak)}; {peak / nnz:.3g} bytes/nz")
    total_measured_time = sum(time.times.values())
    err(f"Total measured time: {total_measured_time:.4g}")
    out_json_path = f"{splitext(memray_bin_path)[0]}.json"
    stats = dict(
        shape=shape,
        nnz=nnz,
        peak_mem=peak,
        times=time.times,
        total_time=total_measured_time,
        bytes_per_nz=peak / nnz,
    )
    json.dump(stats, stdout, indent=2)
    print()
    with open(out_json_path, "w") as f:
        json.dump(stats, f, indent=2)


if __name__ == "__main__":
    cli()
