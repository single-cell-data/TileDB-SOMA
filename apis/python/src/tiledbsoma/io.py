import math
import time
from typing import List, Optional, Sequence, Tuple, Union

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pyarrow as pa
import scipy.sparse as sp
from anndata._core.sparse_dataset import SparseDataset
from somacore.options import PlatformConfig

import tiledbsoma.eta as eta
import tiledbsoma.util as util
import tiledbsoma.util_ann as util_ann
import tiledbsoma.util_tiledb as util_tiledb
from tiledbsoma import (
    Collection,
    DataFrame,
    DenseNDArray,
    Experiment,
    Measurement,
    SparseNDArray,
    logging,
    util_scipy,
)
from tiledbsoma.exception import SOMAError

from .constants import SOMA_JOINID
from .options import SOMATileDBContext, TileDBCreateOptions
from .types import INGEST_MODES, IngestMode, NPNDArray, Path

SparseMatrix = Union[sp.csr_matrix, sp.csc_matrix, SparseDataset]
Matrix = Union[NPNDArray, SparseMatrix]


# ----------------------------------------------------------------
def from_h5ad(
    experiment: Experiment,
    input_path: Path,
    measurement_name: str,
    *,
    context: Optional[SOMATileDBContext] = None,
    platform_config: Optional[PlatformConfig] = None,
    ingest_mode: IngestMode = "write",
) -> None:
    """
    Reads an .h5ad file and writes to a TileDB group structure.

    The "write" ingest_mode (which is the default) writes all data, creating new layers if the soma already exists.

    The "resume" ingest_mode skips data writes if data are within dimension ranges of the existing soma.
    This is useful for continuing after a partial/interrupted previous upload.

    The "schema_only" ingest_mode creates groups and array schema, without writing array data.
    This is useful as a prep-step for parallel append-ingest of multiple H5ADs to a single soma.
    """
    if ingest_mode not in INGEST_MODES:
        raise SOMAError(
            f'expected ingest_mode to be one of {INGEST_MODES}; got "{ingest_mode}"'
        )

    if isinstance(input_path, ad.AnnData):
        raise TypeError("Input path is an AnnData object -- did you want from_anndata?")

    s = util.get_start_stamp()
    logging.log_io(None, f"START  Experiment.from_h5ad {input_path}")

    logging.log_io(None, f"START  READING {input_path}")

    anndata = ad.read_h5ad(input_path, backed="r")

    logging.log_io(None, util.format_elapsed(s, f"FINISH READING {input_path}"))

    from_anndata(
        experiment,
        anndata,
        measurement_name,
        context=context,
        platform_config=platform_config,
        ingest_mode=ingest_mode,
    )

    logging.log_io(
        None, util.format_elapsed(s, f"FINISH Experiment.from_h5ad {input_path}")
    )


# ----------------------------------------------------------------
def from_anndata(
    experiment: Experiment,
    anndata: ad.AnnData,
    measurement_name: str,
    *,
    context: Optional[SOMATileDBContext] = None,
    platform_config: Optional[PlatformConfig] = None,
    ingest_mode: IngestMode = "write",
) -> None:
    """
    Top-level writer method for creating a TileDB group for a ``Experiment`` object.

    The "write" ingest_mode (which is the default) writes all data, creating new layers if the soma already exists.

    The "resume" ingest_mode skips data writes if data are within dimension ranges of the existing soma.
    This is useful for continuing after a partial/interrupted previous upload.

    The "schema_only" ingest_mode creates groups and array schema, without writing array data.
    This is useful as a prep-step for parallel append-ingest of multiple H5ADs to a single soma.
    """
    if ingest_mode not in INGEST_MODES:
        raise SOMAError(
            f'expected ingest_mode to be one of {INGEST_MODES}; got "{ingest_mode}"'
        )

    if not isinstance(anndata, ad.AnnData):
        raise TypeError(
            "Second argument is not an AnnData object -- did you want from_h5ad?"
        )

    # Without _at least_ an index, there is nothing to indicate the dimension indices.
    if anndata.obs.index.empty or anndata.var.index.empty:
        raise NotImplementedError("Empty AnnData.obs or AnnData.var unsupported.")

    s = util.get_start_stamp()
    logging.log_io(None, "START  DECATEGORICALIZING")

    anndata.obs_names_make_unique()
    anndata.var_names_make_unique()

    logging.log_io(None, util.format_elapsed(s, "FINISH DECATEGORICALIZING"))

    s = util.get_start_stamp()
    logging.log_io(None, f"START  WRITING {experiment.uri}")

    # Must be done first, to create the parent directory.
    _check_create_experiment(experiment, ingest_mode)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # OBS
    obs = DataFrame(uri=util.uri_joinpath(experiment.uri, "obs"))
    _write_dataframe(
        obs,
        util_ann._decategoricalize_obs_or_var(anndata.obs),
        id_column_name="obs_id",
        platform_config=platform_config,
        ingest_mode=ingest_mode,
    )
    experiment.set("obs", obs)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # MS
    ms = Collection(uri=util.uri_joinpath(experiment.uri, "ms"))

    experiment.set("ms", _check_create_collection(ms, ingest_mode))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # MS/meas
    measurement = _check_create_measurement(
        Measurement(uri=f"{experiment.ms.uri}/{measurement_name}", context=context),
        ingest_mode,
    )
    experiment.ms.set(measurement_name, measurement)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # MS/meas/VAR
    var = DataFrame(uri=util.uri_joinpath(measurement.uri, "var"))
    _write_dataframe(
        var,
        util_ann._decategoricalize_obs_or_var(anndata.var),
        id_column_name="var_id",
        platform_config=platform_config,
        ingest_mode=ingest_mode,
    )
    measurement["var"] = var

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # MS/meas/X/DATA

    x = Collection(uri=util.uri_joinpath(measurement.uri, "X"))
    measurement["X"] = _check_create_collection(x, ingest_mode)

    # Since we did `anndata = ad.read_h5ad(path_to_h5ad, "r")` with the "r":
    # * If we do `anndata.X[:]` we're loading all of a CSR/CSC/etc into memory.
    # * If we do `anndata.X` we're getting a pageable object which can be loaded
    #   chunkwise into memory.
    # Using the latter allows us to ingest larger .h5ad files without OOMing.
    cls = (
        DenseNDArray
        if isinstance(anndata.X, (np.ndarray, h5py.Dataset))
        else SparseNDArray
    )
    data = cls(uri=util.uri_joinpath(measurement.X.uri, "data"), context=context)
    create_from_matrix(data, anndata.X, platform_config, ingest_mode)
    measurement.X.set("data", data)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # MS/meas/OBSM,VARM,OBSP,VARP
    if len(anndata.obsm.keys()) > 0:  # do not create an empty collection
        measurement["obsm"] = _check_create_collection(
            Collection(uri=util.uri_joinpath(measurement.uri, "obsm")),
            ingest_mode,
        )
        for key in anndata.obsm.keys():
            arr = DenseNDArray(
                uri=util.uri_joinpath(measurement.obsm.uri, key),
                context=context,
            )
            create_from_matrix(
                arr,
                util_tiledb.to_tiledb_supported_array_type(anndata.obsm[key]),
                platform_config,
                ingest_mode,
            )
            measurement.obsm.set(key, arr)

    if len(anndata.varm.keys()) > 0:  # do not create an empty collection
        measurement["varm"] = _check_create_collection(
            Collection(uri=util.uri_joinpath(measurement.uri, "varm")),
            ingest_mode,
        )
        for key in anndata.varm.keys():
            darr = DenseNDArray(
                uri=util.uri_joinpath(measurement.varm.uri, key),
                context=context,
            )
            create_from_matrix(
                darr,
                util_tiledb.to_tiledb_supported_array_type(anndata.varm[key]),
                platform_config,
                ingest_mode,
            )
            measurement.varm.set(key, darr)

    if len(anndata.obsp.keys()) > 0:  # do not create an empty collection
        measurement["obsp"] = _check_create_collection(
            Collection(uri=util.uri_joinpath(measurement.uri, "obsp")),
            ingest_mode,
        )
        for key in anndata.obsp.keys():
            sarr = SparseNDArray(
                uri=util.uri_joinpath(measurement.obsp.uri, key),
                context=context,
            )
            create_from_matrix(
                sarr,
                util_tiledb.to_tiledb_supported_array_type(anndata.obsp[key]),
                platform_config,
                ingest_mode,
            )
            measurement.obsp.set(key, sarr)

    if len(anndata.varp.keys()) > 0:  # do not create an empty collection
        measurement["varp"] = _check_create_collection(
            Collection(uri=util.uri_joinpath(measurement.uri, "varp")),
            ingest_mode,
        )
        for key in anndata.varp.keys():
            sarr = SparseNDArray(
                uri=util.uri_joinpath(measurement.varp.uri, key),
                context=context,
            )
            create_from_matrix(
                sarr,
                util_tiledb.to_tiledb_supported_array_type(anndata.varp[key]),
                platform_config,
                ingest_mode,
            )
            measurement.varp.set(key, sarr)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # MS/RAW
    if anndata.raw is not None:
        raw_measurement = Measurement(
            uri=util.uri_joinpath(experiment.ms.uri, "raw"),
            context=context,
        )
        experiment.ms.set(
            "raw", _check_create_measurement(raw_measurement, ingest_mode)
        )

        var = DataFrame(uri=util.uri_joinpath(raw_measurement.uri, "var"))
        _write_dataframe(
            var,
            util_ann._decategoricalize_obs_or_var(anndata.raw.var),
            id_column_name="var_id",
            platform_config=platform_config,
            ingest_mode=ingest_mode,
        )
        raw_measurement.set("var", var)

        raw_measurement["X"] = _check_create_collection(
            Collection(uri=util.uri_joinpath(raw_measurement.uri, "X")),
            ingest_mode,
        )

        rawXdata = SparseNDArray(
            uri=util.uri_joinpath(raw_measurement.X.uri, "data"),
            context=context,
        )
        create_from_matrix(rawXdata, anndata.raw.X, platform_config, ingest_mode)
        raw_measurement.X.set("data", rawXdata)

    logging.log_io(
        f"Wrote   {experiment.uri}",
        util.format_elapsed(s, f"FINISH WRITING {experiment.uri}"),
    )


def _check_create(
    thing: Union[Experiment, Measurement, Collection],
    ingest_mode: str,
) -> Union[Experiment, Measurement, Collection]:
    if ingest_mode == "resume":
        if not thing.exists():
            thing.create()
        # else fine
    else:
        if thing.exists():
            raise SOMAError(f"{thing.uri} already exists")
        else:
            thing.create()
    return thing


# Split out from _check_create since its union-return gives the type-checker fits at callsites.
def _check_create_collection(thing: Collection, ingest_mode: str) -> Collection:
    retval = _check_create(thing, ingest_mode)
    if not isinstance(retval, Collection):
        raise SOMAError(
            f"internal error: expected object of type Collection; got {type(retval)}"
        )

    return retval


def _check_create_experiment(thing: Experiment, ingest_mode: str) -> Experiment:
    retval = _check_create(thing, ingest_mode)
    if not isinstance(retval, Experiment):
        raise SOMAError(
            f"internal error: expected object of type Experiment; got {type(retval)}"
        )
    return retval


def _check_create_measurement(thing: Measurement, ingest_mode: str) -> Measurement:
    retval = _check_create(thing, ingest_mode)
    if not isinstance(retval, Measurement):
        raise SOMAError(
            f"internal error: expected object of type Measurement; got {type(retval)}"
        )
    return retval


def _write_dataframe(
    soma_df: DataFrame,
    df: pd.DataFrame,
    id_column_name: Optional[str],
    platform_config: Optional[PlatformConfig] = None,
    ingest_mode: IngestMode = "write",
) -> None:
    s = util.get_start_stamp()
    logging.log_io(None, f"START  WRITING {soma_df.uri}")

    df[SOMA_JOINID] = np.asarray(range(len(df)), dtype=np.int64)

    df.reset_index(inplace=True)
    if id_column_name is not None:
        df.rename(columns={"index": id_column_name}, inplace=True)
    df.set_index(SOMA_JOINID, inplace=True)

    # Categoricals are not yet well supported, so we must flatten
    for k in df:
        if df[k].dtype == "category":
            df[k] = df[k].astype(df[k].cat.categories.dtype)
    arrow_table = pa.Table.from_pandas(df)

    if soma_df.exists():
        if ingest_mode == "resume":
            # This lets us check for already-ingested dataframes, when in resume-ingest mode.
            with soma_df._tiledb_open() as A:
                storage_ned = A.nonempty_domain()
            dim_range = ((int(df.index.min()), int(df.index.max())),)
            if _chunk_is_contained_in(dim_range, storage_ned):
                logging.log_io(
                    f"Skipped {soma_df.uri}",
                    util.format_elapsed(s, f"SKIPPED {soma_df.uri}"),
                )
                return
        else:
            raise SOMAError(f"{soma_df.uri} already exists")
    else:
        soma_df.create(arrow_table.schema, platform_config=platform_config)

    if ingest_mode == "schema_only":
        logging.log_io(
            f"Wrote schema {soma_df.uri}",
            util.format_elapsed(s, f"FINISH WRITING SCHEMA {soma_df.uri}"),
        )
        return

    soma_df.write(arrow_table)
    logging.log_io(
        f"Wrote   {soma_df.uri}",
        util.format_elapsed(s, f"FINISH WRITING {soma_df.uri}"),
    )


def create_from_matrix(
    soma_ndarray: Union[DenseNDArray, SparseNDArray],
    matrix: Union[Matrix, h5py.Dataset],
    platform_config: Optional[PlatformConfig] = None,
    ingest_mode: IngestMode = "write",
) -> None:
    """
    Create and populate the ``soma_matrix`` from the contents of ``matrix``.
    """
    # SparseDataset has no ndim but it has a shape
    if len(matrix.shape) != 2:
        raise ValueError(f"expected matrix.shape == 2; got {matrix.shape}")
    acceptables = ("SOMADenseNDArray", "SOMASparseNDArray")
    if soma_ndarray.soma_type not in acceptables:
        raise ValueError(
            f'internal error: expected array type to be one of {acceptables}; got "{soma_ndarray.soma_type}"'
        )

    s = util.get_start_stamp()
    logging.log_io(None, f"START  WRITING {soma_ndarray.uri}")

    if ingest_mode != "resume" or not soma_ndarray.exists():
        if soma_ndarray.exists():
            raise SOMAError(f"{soma_ndarray.uri} already exists")
        soma_ndarray.create(
            type=pa.from_numpy_dtype(matrix.dtype),
            shape=matrix.shape,
            platform_config=platform_config,
        )

    if ingest_mode == "schema_only":
        logging.log_io(
            f"Wrote schema {soma_ndarray.uri}",
            util.format_elapsed(s, f"FINISH WRITING SCHEMA {soma_ndarray.uri}"),
        )
        return

    logging.log_io(
        f"Writing {soma_ndarray.uri}",
        util.format_elapsed(s, f"START  WRITING {soma_ndarray.uri}"),
    )

    if isinstance(soma_ndarray, DenseNDArray):
        _write_matrix_to_denseNDArray(
            soma_ndarray,
            matrix,
            tiledb_create_options=TileDBCreateOptions.from_platform_config(
                platform_config
            ),
            ingest_mode=ingest_mode,
        )
    else:  # SOMASparseNDArray
        _write_matrix_to_sparseNDArray(
            soma_ndarray,
            matrix,
            tiledb_create_options=TileDBCreateOptions.from_platform_config(
                platform_config
            ),
            ingest_mode=ingest_mode,
        )

    logging.log_io(
        f"Wrote   {soma_ndarray.uri}",
        util.format_elapsed(s, f"FINISH WRITING {soma_ndarray.uri}"),
    )


def _write_matrix_to_denseNDArray(
    soma_ndarray: DenseNDArray,
    matrix: Union[Matrix, h5py.Dataset],
    tiledb_create_options: TileDBCreateOptions,
    ingest_mode: IngestMode,
) -> None:
    """Write a matrix to an empty DenseNDArray"""

    # There is a chunk-by-chunk already-done check for resume mode, below.
    # This full-matrix-level check here might seem redundant, but in fact it's important:
    # * By checking input bounds against storage NED here, we can see if the entire matrix
    #   was already ingested and avoid even loading chunks;
    # * By checking chunkwise we can catch the case where a matrix was already *partly*
    #   ingested.
    # * Of course, this also helps us catch already-completed writes in the non-chunked case.
    storage_ned = None
    if ingest_mode == "resume" and soma_ndarray.exists():
        # This lets us check for already-ingested chunks, when in resume-ingest mode.
        with soma_ndarray._tiledb_open() as A:
            storage_ned = A.nonempty_domain()
            matrix_bounds = [
                (0, int(n - 1)) for n in matrix.shape
            ]  # Cast for lint in case np.int64
            logging.log_io(
                None,
                f"Input bounds {tuple(matrix_bounds)} storage non-empty domain {storage_ned}",
            )
            if _chunk_is_contained_in(matrix_bounds, storage_ned):
                logging.log_io(
                    f"Skipped {soma_ndarray.uri}", f"SKIPPED WRITING {soma_ndarray.uri}"
                )
                return

    # Write all at once?
    if not tiledb_create_options.write_X_chunked():
        if not isinstance(matrix, np.ndarray):
            matrix = matrix.toarray()
        soma_ndarray.write((slice(None),), pa.Tensor.from_numpy(matrix))
        return

    # OR, write in chunks
    eta_tracker = eta.Tracker()
    nrow, ncol = matrix.shape
    i = 0
    # Number of rows to chunk by. Dense writes, so this is a constant.
    chunk_size = int(math.ceil(tiledb_create_options.goal_chunk_nnz() / ncol))
    while i < nrow:
        t1 = time.time()
        i2 = i + chunk_size

        # Print doubly-inclusive lo..hi like 0..17 and 18..31.
        chunk_percent = min(100, 100 * (i2 - 1) / nrow)
        logging.log_io(
            None,
            "START  chunk rows %d..%d of %d (%.3f%%)"
            % (i, i2 - 1, nrow, chunk_percent),
        )

        chunk = matrix[i:i2, :]

        if ingest_mode == "resume" and storage_ned is not None:
            chunk_bounds = matrix_bounds
            chunk_bounds[0] = (int(i), int(i2 - 1))  # Cast for lint in case np.int64
            if _chunk_is_contained_in_axis(chunk_bounds, storage_ned, 0):
                # Print doubly inclusive lo..hi like 0..17 and 18..31.
                logging.log_io(
                    "... %7.3f%% done" % chunk_percent,
                    "SKIP   chunk rows %d..%d of %d (%.3f%%)"
                    % (i, i2 - 1, nrow, chunk_percent),
                )
                i = i2
                continue

        if isinstance(chunk, np.ndarray):
            tensor = pa.Tensor.from_numpy(chunk)
        else:
            tensor = pa.Tensor.from_numpy(chunk.toarray())
        soma_ndarray.write((slice(i, i2), slice(None)), tensor)

        t2 = time.time()
        chunk_seconds = t2 - t1
        eta_seconds = eta_tracker.ingest_and_predict(chunk_percent, chunk_seconds)

        if chunk_percent < 100:
            logging.log_io(
                "... %7.3f%% done, ETA %s" % (chunk_percent, eta_seconds),
                "FINISH chunk in %.3f seconds, %7.3f%% done, ETA %s"
                % (chunk_seconds, chunk_percent, eta_seconds),
            )

        i = i2

    return


def _find_sparse_chunk_size(
    matrix: SparseMatrix, start_index: int, axis: int, goal_chunk_nnz: int
) -> int:
    """
    Given a sparse matrix and a start index, return a step size, on the stride axis, which will
    achieve the cummulative nnz desired.

    :param matrix: The input scipy.sparse matrix.
    :param start_index: the index at which to start a chunk.
    :param axis: the stride axis, across which to find a chunk.
    :param goal_chunk_nnz: Desired number of non-zero array entries for the chunk.
    """
    chunk_size = 1
    sum_nnz = 0
    coords: List[Union[slice, int]] = [slice(None), slice(None)]

    # Empirically we find:
    # * If the input matrix is sp.csr_matrix or sp.csc_matrix then getting all these nnz values is
    #   quick.
    # * If the input matrix is anndata._core.sparse_dataset.SparseDataset -- which happens with
    #   out-of-core anndata reads -- then getting all these nnz values is prohibitively expensive.
    # * It turns out that getting a sample is quite sufficient. We do this regardless of whether
    #   the matrix is anndata._core.sparse_dataset.SparseDataset or not.
    # * The max_rows is manually defined after running experiments with 60GB .h5ad files.
    count = 0
    max_rows = 100

    for index in range(start_index, matrix.shape[axis]):
        count += 1
        coords[axis] = index
        sum_nnz += matrix[tuple(coords)].nnz
        if sum_nnz > goal_chunk_nnz:
            break
        if count > max_rows:
            break
        chunk_size += 1

    if sum_nnz > goal_chunk_nnz:
        return chunk_size

    # Solve the equation:
    #
    # sum_nnz              count
    # -------          =  -------
    # goal_chunk_nnz       result
    chunk_size = int(count * goal_chunk_nnz / sum_nnz)
    if chunk_size < 1:
        chunk_size = 1
    return chunk_size


def _write_matrix_to_sparseNDArray(
    soma_ndarray: SparseNDArray,
    matrix: Matrix,
    tiledb_create_options: TileDBCreateOptions,
    ingest_mode: IngestMode,
) -> None:
    """Write a matrix to an empty DenseNDArray"""

    def _coo_to_table(mat_coo: sp.coo_matrix, axis: int = 0, base: int = 0) -> pa.Table:
        pydict = {
            "soma_data": mat_coo.data,
            "soma_dim_0": mat_coo.row + base if base > 0 and axis == 0 else mat_coo.row,
            "soma_dim_1": mat_coo.col + base if base > 0 and axis == 1 else mat_coo.col,
        }
        return pa.Table.from_pydict(pydict)

    # There is a chunk-by-chunk already-done check for resume mode, below.
    # This full-matrix-level check here might seem redundant, but in fact it's important:
    # * By checking input bounds against storage NED here, we can see if the entire matrix
    #   was already ingested and avoid even loading chunks;
    # * By checking chunkwise we can catch the case where a matrix was already *partly*
    #   ingested.
    # * Of course, this also helps us catch already-completed writes in the non-chunked case.
    storage_ned = None
    if ingest_mode == "resume" and soma_ndarray.exists():
        # This lets us check for already-ingested chunks, when in resume-ingest mode.
        with soma_ndarray._tiledb_open() as A:
            storage_ned = A.nonempty_domain()
            matrix_bounds = [
                (0, int(n - 1)) for n in matrix.shape
            ]  # Cast for lint in case np.int64
            logging.log_io(
                None,
                f"Input bounds {tuple(matrix_bounds)} storage non-empty domain {storage_ned}",
            )
            if _chunk_is_contained_in(matrix_bounds, storage_ned):
                logging.log_io(
                    f"Skipped {soma_ndarray.uri}", f"SKIPPED WRITING {soma_ndarray.uri}"
                )
                return

    # Write all at once?
    if not tiledb_create_options.write_X_chunked():
        soma_ndarray.write(_coo_to_table(sp.coo_matrix(matrix)))
        return

    # Or, write in chunks, striding across the most efficient slice axis

    stride_axis = 0
    if sp.isspmatrix_csc(matrix):
        # E.g. if we used anndata.X[:]
        stride_axis = 1
    if isinstance(matrix, SparseDataset) and matrix.format_str == "csc":
        # E.g. if we used anndata.X without the [:]
        stride_axis = 1

    dim_max_size = matrix.shape[stride_axis]

    eta_tracker = eta.Tracker()
    goal_chunk_nnz = tiledb_create_options.goal_chunk_nnz()

    coords = [slice(None), slice(None)]
    i = 0
    while i < dim_max_size:
        t1 = time.time()

        # Chunk size on the stride axis
        if isinstance(matrix, np.ndarray):
            chunk_size = int(math.ceil(goal_chunk_nnz / matrix.shape[stride_axis]))
        else:
            chunk_size = _find_sparse_chunk_size(matrix, i, stride_axis, goal_chunk_nnz)

        i2 = i + chunk_size

        coords[stride_axis] = slice(i, i2)
        chunk_coo = sp.coo_matrix(matrix[tuple(coords)])

        chunk_percent = min(100, 100 * (i2 - 1) / dim_max_size)

        if ingest_mode == "resume" and storage_ned is not None:
            chunk_bounds = matrix_bounds
            chunk_bounds[stride_axis] = (
                int(i),
                int(i2 - 1),
            )  # Cast for lint in case np.int64
            if _chunk_is_contained_in_axis(chunk_bounds, storage_ned, stride_axis):
                # Print doubly inclusive lo..hi like 0..17 and 18..31.
                logging.log_io(
                    "... %7.3f%% done" % chunk_percent,
                    "SKIP   chunk rows %d..%d of %d (%.3f%%), nnz=%d"
                    % (i, i2 - 1, dim_max_size, chunk_percent, chunk_coo.nnz),
                )
                i = i2
                continue

        # Print doubly inclusive lo..hi like 0..17 and 18..31.
        logging.log_io(
            None,
            "START  chunk rows %d..%d of %d (%.3f%%), nnz=%d"
            % (i, i2 - 1, dim_max_size, chunk_percent, chunk_coo.nnz),
        )

        soma_ndarray.write(_coo_to_table(chunk_coo, stride_axis, i))

        t2 = time.time()
        chunk_seconds = t2 - t1
        eta_seconds = eta_tracker.ingest_and_predict(chunk_percent, chunk_seconds)

        if chunk_percent < 100:
            logging.log_io(
                "... %7.3f%% done, ETA %s" % (chunk_percent, eta_seconds),
                "FINISH chunk in %.3f seconds, %7.3f%% done, ETA %s"
                % (chunk_seconds, chunk_percent, eta_seconds),
            )

        i = i2


def _chunk_is_contained_in(
    chunk_bounds: Sequence[Tuple[int, int]],
    storage_nonempty_domain: Optional[Sequence[Tuple[Optional[int], Optional[int]]]],
) -> bool:
    """
    Determines if a dim range is included within the array's non-empty domain.  Ranges are inclusive
    on both endpoints.  This is a helper for resume-ingest mode.

    We say "bounds" not "MBR" with the "M" for minimum: a sparse matrix might not _have_ any
    elements for some initial/final rows or columns. Suppose an input array has shape 100 x 200, so
    bounds `((0, 99), (0, 199))` -- and also suppose there are no matrix elements for column 1.
    Also suppose the matrix has already been written to TileDB-SOMA storage. The TileDB non-empty
    domain _is_ tight -- it'd say `((0, 99), (3, 197))` for example.  When we come back for a
    resume-mode ingest, we'd see the input bounds aren't contained within the storage non-empty
    domain, and erroneously declare that the data need to be rewritten.

    This is why we take the stride axis as an argument. In resume mode, it's our contract with the
    user that they declare they are retrying the exact same input file -- and we do our best to
    fulfill their ask by checking the dimension being strided on.
    """
    if storage_nonempty_domain is None:
        return False

    if len(chunk_bounds) != len(storage_nonempty_domain):
        raise SOMAError(
            f"internal error: ingest data ndim {len(chunk_bounds)} != storage ndim {len(storage_nonempty_domain)}"
        )
    for i in range(len(chunk_bounds)):
        if not _chunk_is_contained_in_axis(chunk_bounds, storage_nonempty_domain, i):
            return False
    return True


def _chunk_is_contained_in_axis(
    chunk_bounds: Sequence[Tuple[int, int]],
    storage_nonempty_domain: Sequence[Tuple[Optional[int], Optional[int]]],
    stride_axis: int,
) -> bool:
    """
    Helper function for ``_chunk_is_contained_in``.
    """
    storage_lo, storage_hi = storage_nonempty_domain[stride_axis]
    if storage_lo is None or storage_hi is None:
        # E.g. an array has had its schema created but no data written yet
        return False

    chunk_lo, chunk_hi = chunk_bounds[stride_axis]
    if chunk_lo < storage_lo or chunk_lo > storage_hi:
        return False
    if chunk_hi < storage_lo or chunk_hi > storage_hi:
        return False

    return True


# ----------------------------------------------------------------
def to_h5ad(
    experiment: Experiment,
    h5ad_path: Path,
    measurement_name: str,
    X_layer_name: str = "data",
) -> None:
    """
    Converts the experiment group to anndata format and writes it to the specified .h5ad file.
    """
    s = util.get_start_stamp()
    logging.log_io(None, f"START  Experiment.to_h5ad -> {h5ad_path}")

    anndata = to_anndata(
        experiment, measurement_name=measurement_name, X_layer_name=X_layer_name
    )

    s2 = util.get_start_stamp()
    logging.log_io(None, f"START  write {h5ad_path}")

    anndata.write_h5ad(h5ad_path)

    logging.log_io(None, util.format_elapsed(s2, f"FINISH write {h5ad_path}"))

    logging.log_io(
        None, util.format_elapsed(s, f"FINISH Experiment.to_h5ad -> {h5ad_path}")
    )


# ----------------------------------------------------------------
def to_anndata(
    experiment: Experiment, measurement_name: str, X_layer_name: str = "data"
) -> ad.AnnData:
    """
    Converts the experiment group to anndata. Choice of matrix formats is following what we often see in input .h5ad files:

    * X as ``scipy.sparse.csr_matrix``
    * obs,var as ``pandas.dataframe``
    * obsm,varm arrays as ``numpy.ndarray``
    * obsp,varp arrays as ``scipy.sparse.csr_matrix``
    """

    s = util.get_start_stamp()
    logging.log_io(None, "START  Experiment.to_anndata")

    measurement: Measurement = experiment.ms[measurement_name]

    obs_df = experiment.obs.read().concat().to_pandas()
    obs_df.drop([SOMA_JOINID], axis=1, inplace=True)
    obs_df.set_index("obs_id", inplace=True)

    var_df = measurement.var.read().concat().to_pandas()
    var_df.drop([SOMA_JOINID], axis=1, inplace=True)
    var_df.set_index("var_id", inplace=True)

    nobs = len(obs_df.index)
    nvar = len(var_df.index)

    if X_layer_name not in measurement.X:
        raise SOMAError(
            f"X_layer_name {X_layer_name} not found in data: {measurement.X.keys()}"
        )
    X_data = measurement.X[X_layer_name]
    X_csr = None
    X_dtype = None  # some datasets have no X
    if isinstance(X_data, DenseNDArray):
        X_ndarray = X_data.read((slice(None), slice(None))).to_numpy()
        X_dtype = X_ndarray.dtype
    elif isinstance(X_data, SparseNDArray):
        X_mat = X_data.read().tables().concat().to_pandas()  # TODO: CSR/CSC options ...
        X_csr = util_scipy.csr_from_tiledb_df(X_mat, nobs, nvar)
        X_dtype = X_csr.dtype
    else:
        raise TypeError(f"Unexpected NDArray type {type(X_data)}")

    obsm = {}
    if "obsm" in measurement and measurement.obsm.exists():
        for key in measurement.obsm.keys():
            shape = measurement.obsm[key].shape
            if len(shape) != 2:
                raise ValueError(f"expected shape == 2; got {shape}")
            matrix = measurement.obsm[key].read((slice(None),) * len(shape)).to_numpy()
            # The spelling `sp.csr_array` is more idiomatic but doesn't exist until Python 3.8
            obsm[key] = sp.csr_matrix(matrix)

    varm = {}
    if "varm" in measurement and measurement.varm.exists():
        for key in measurement.varm.keys():
            shape = measurement.varm[key].shape
            if len(shape) != 2:
                raise ValueError(f"expected shape == 2; got {shape}")
            matrix = measurement.varm[key].read((slice(None),) * len(shape)).to_numpy()
            # The spelling `sp.csr_array` is more idiomatic but doesn't exist until Python 3.8
            varm[key] = sp.csr_matrix(matrix)

    obsp = {}
    if "obsp" in measurement and measurement.obsp.exists():
        for key in measurement.obsp.keys():
            matrix = measurement.obsp[key].read().tables().concat().to_pandas()
            obsp[key] = util_scipy.csr_from_tiledb_df(matrix, nobs, nobs)

    varp = {}
    if "varp" in measurement and measurement.varp.exists():
        for key in measurement.varp.keys():
            matrix = measurement.varp[key].read().tables().concat().to_pandas()
            varp[key] = util_scipy.csr_from_tiledb_df(matrix, nvar, nvar)

    anndata = ad.AnnData(
        X=X_csr if X_csr is not None else X_ndarray,
        obs=obs_df,
        var=var_df,
        obsm=obsm,
        varm=varm,
        obsp=obsp,
        varp=varp,
        dtype=X_dtype,
    )

    logging.log_io(None, util.format_elapsed(s, "FINISH Experiment.to_anndata"))

    return anndata
