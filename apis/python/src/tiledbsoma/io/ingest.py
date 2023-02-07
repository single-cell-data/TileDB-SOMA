import math
import time
from typing import (
    Any,
    List,
    Optional,
    Sequence,
    Tuple,
    Type,
    TypeVar,
    Union,
    cast,
    overload,
)

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pyarrow as pa
import scipy.sparse as sp
from anndata._core.sparse_dataset import SparseDataset
from somacore.options import PlatformConfig

from .. import (
    Collection,
    DataFrame,
    DenseNDArray,
    Experiment,
    Measurement,
    SparseNDArray,
    eta,
    factory,
    logging,
    util,
)
from ..common_nd_array import NDArray
from ..constants import SOMA_JOINID
from ..exception import DoesNotExistError, SOMAError
from ..funcs import typeguard_ignore
from ..options import SOMATileDBContext
from ..options.tiledb_create_options import TileDBCreateOptions
from ..tdb_handles import RawHandle
from ..tiledb_array import TileDBArray
from ..tiledb_object import TileDBObject
from ..types import INGEST_MODES, IngestMode, NPNDArray, Path
from . import conversions

SparseMatrix = Union[sp.csr_matrix, sp.csc_matrix, SparseDataset]
Matrix = Union[NPNDArray, SparseMatrix]
_NDArr = TypeVar("_NDArr", bound=NDArray)
_TDBO = TypeVar("_TDBO", bound=TileDBObject[RawHandle])


# ----------------------------------------------------------------
def from_h5ad(
    experiment_uri: str,
    input_path: Path,
    measurement_name: str,
    *,
    context: Optional[SOMATileDBContext] = None,
    platform_config: Optional[PlatformConfig] = None,
    ingest_mode: IngestMode = "write",
    use_relative_uri: Optional[bool] = None,
) -> Experiment:
    """
    Reads an .h5ad file and writes to a TileDB group structure.

    Returns an experiment opened for writing.

    The "write" ingest_mode (which is the default) writes all data, creating new layers if the soma already exists.

    The "resume" ingest_mode skips data writes if data are within dimension ranges of the existing soma.
    This is useful for continuing after a partial/interrupted previous upload.

    The "schema_only" ingest_mode creates groups and array schema, without writing array data.
    This is useful as a prep-step for parallel append-ingest of multiple H5ADs to a single soma.

    [lifecycle: experimental]
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

    exp = from_anndata(
        experiment_uri,
        anndata,
        measurement_name,
        context=context,
        platform_config=platform_config,
        ingest_mode=ingest_mode,
        use_relative_uri=use_relative_uri,
    )

    logging.log_io(
        None, util.format_elapsed(s, f"FINISH Experiment.from_h5ad {input_path}")
    )
    return exp


# ----------------------------------------------------------------
def from_anndata(
    experiment_uri: str,
    anndata: ad.AnnData,
    measurement_name: str,
    *,
    context: Optional[SOMATileDBContext] = None,
    platform_config: Optional[PlatformConfig] = None,
    ingest_mode: IngestMode = "write",
    use_relative_uri: Optional[bool] = None,
) -> Experiment:
    """
    Top-level writer method for creating a TileDB group for a ``Experiment`` object.

    Returns an Experiment opened for writing.

    The "write" ingest_mode (which is the default) writes all data, creating new layers if the soma already exists.

    The "resume" ingest_mode skips data writes if data are within dimension ranges of the existing soma.
    This is useful for continuing after a partial/interrupted previous upload.

    The "schema_only" ingest_mode creates groups and array schema, without writing array data.
    This is useful as a prep-step for parallel append-ingest of multiple H5ADs to a single soma.

    [lifecycle: experimental]
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
    logging.log_io(None, f"START  WRITING {experiment_uri}")

    # Must be done first, to create the parent directory.
    experiment = _create_or_open_coll(Experiment, experiment_uri, ingest_mode)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # OBS
    df_uri = util.uri_joinpath(experiment.uri, "obs")
    with _write_dataframe(
        df_uri,
        conversions.decategoricalize_obs_or_var(anndata.obs),
        id_column_name="obs_id",
        platform_config=platform_config,
        ingest_mode=ingest_mode,
    ) as obs:
        experiment.set("obs", obs, use_relative_uri=use_relative_uri)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # MS
    with _create_or_open_coll(
        Collection[Measurement], util.uri_joinpath(experiment.uri, "ms"), ingest_mode
    ) as ms:
        experiment.set("ms", ms, use_relative_uri=use_relative_uri)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # MS/meas
        with _create_or_open_coll(
            Measurement, f"{experiment.ms.uri}/{measurement_name}", ingest_mode
        ) as measurement:
            ms.set(measurement_name, measurement, use_relative_uri=use_relative_uri)

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # MS/meas/VAR
            with _write_dataframe(
                util.uri_joinpath(measurement.uri, "var"),
                conversions.decategoricalize_obs_or_var(anndata.var),
                id_column_name="var_id",
                platform_config=platform_config,
                ingest_mode=ingest_mode,
            ) as var:
                measurement.set("var", var, use_relative_uri=use_relative_uri)

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # MS/meas/X/DATA

            with _create_or_open_coll(
                Collection, util.uri_joinpath(measurement.uri, "X"), ingest_mode
            ) as x:
                measurement.set("X", x, use_relative_uri=use_relative_uri)

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
                with create_from_matrix(
                    cls,
                    util.uri_joinpath(measurement.X.uri, "data"),
                    anndata.X,
                    platform_config,
                    ingest_mode,
                ) as data:
                    x.set("data", data, use_relative_uri=use_relative_uri)

                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                # MS/meas/OBSM,VARM,OBSP,VARP
                if len(anndata.obsm.keys()) > 0:  # do not create an empty collection
                    with _create_or_open_coll(
                        Collection,
                        util.uri_joinpath(measurement.uri, "obsm"),
                        ingest_mode,
                    ) as obsm:
                        measurement.set("obsm", obsm, use_relative_uri=use_relative_uri)
                        for key in anndata.obsm.keys():
                            with create_from_matrix(
                                DenseNDArray,
                                util.uri_joinpath(measurement.obsm.uri, key),
                                conversions.to_tiledb_supported_array_type(
                                    anndata.obsm[key]
                                ),
                                platform_config,
                                ingest_mode,
                            ) as arr:
                                obsm.set(key, arr, use_relative_uri=use_relative_uri)
                            arr.close()
                    measurement.obsm.close()

                if len(anndata.varm.keys()) > 0:  # do not create an empty collection
                    with _create_or_open_coll(
                        Collection,
                        util.uri_joinpath(measurement.uri, "varm"),
                        ingest_mode,
                    ) as varm:
                        measurement.set("varm", varm, use_relative_uri=use_relative_uri)
                        for key in anndata.varm.keys():
                            with create_from_matrix(
                                DenseNDArray,
                                util.uri_joinpath(measurement.varm.uri, key),
                                conversions.to_tiledb_supported_array_type(
                                    anndata.varm[key]
                                ),
                                platform_config,
                                ingest_mode,
                            ) as darr:
                                varm.set(
                                    key,
                                    darr,
                                    use_relative_uri=use_relative_uri,
                                )

                if len(anndata.obsp.keys()) > 0:  # do not create an empty collection
                    with _create_or_open_coll(
                        Collection,
                        util.uri_joinpath(measurement.uri, "obsp"),
                        ingest_mode,
                    ) as obsp:
                        measurement.set("obsp", obsp, use_relative_uri=use_relative_uri)
                        for key in anndata.obsp.keys():
                            with create_from_matrix(
                                SparseNDArray,
                                util.uri_joinpath(measurement.obsp.uri, key),
                                conversions.to_tiledb_supported_array_type(
                                    anndata.obsp[key]
                                ),
                                platform_config,
                                ingest_mode,
                            ) as sarr:
                                obsp.set(
                                    key,
                                    sarr,
                                    use_relative_uri=use_relative_uri,
                                )

                if len(anndata.varp.keys()) > 0:  # do not create an empty collection
                    with _create_or_open_coll(
                        Collection,
                        util.uri_joinpath(measurement.uri, "varp"),
                        ingest_mode,
                    ) as varp:
                        measurement.set("varp", varp, use_relative_uri=use_relative_uri)
                        for key in anndata.varp.keys():
                            with create_from_matrix(
                                SparseNDArray,
                                util.uri_joinpath(measurement.varp.uri, key),
                                conversions.to_tiledb_supported_array_type(
                                    anndata.varp[key]
                                ),
                                platform_config,
                                ingest_mode,
                            ) as sarr:
                                varp.set(
                                    key,
                                    sarr,
                                    use_relative_uri=use_relative_uri,
                                )

                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                # MS/RAW
                if anndata.raw is not None:
                    with _create_or_open_coll(
                        Measurement,
                        util.uri_joinpath(experiment.ms.uri, "raw"),
                        ingest_mode,
                    ) as raw_measurement:
                        ms.set(
                            "raw",
                            raw_measurement,
                            use_relative_uri=use_relative_uri,
                        )

                        with _write_dataframe(
                            util.uri_joinpath(raw_measurement.uri, "var"),
                            conversions.decategoricalize_obs_or_var(anndata.raw.var),
                            id_column_name="var_id",
                            platform_config=platform_config,
                            ingest_mode=ingest_mode,
                        ) as var:
                            raw_measurement.set(
                                "var", var, use_relative_uri=use_relative_uri
                            )

                        with _create_or_open_coll(
                            Collection,
                            util.uri_joinpath(raw_measurement.uri, "X"),
                            ingest_mode,
                        ) as rm_x:
                            raw_measurement.set(
                                "X", rm_x, use_relative_uri=use_relative_uri
                            )

                            with create_from_matrix(
                                SparseNDArray,
                                util.uri_joinpath(raw_measurement.X.uri, "data"),
                                anndata.raw.X,
                                platform_config,
                                ingest_mode,
                            ) as rm_x_data:
                                rm_x.set(
                                    "data",
                                    rm_x_data,
                                    use_relative_uri=use_relative_uri,
                                )

    logging.log_io(
        f"Wrote   {experiment.uri}",
        util.format_elapsed(s, f"FINISH WRITING {experiment.uri}"),
    )
    return experiment


@overload
def _create_or_open_coll(
    cls: Type[Experiment], uri: str, ingest_mode: str
) -> Experiment:
    ...


@overload
def _create_or_open_coll(
    cls: Type[Measurement], uri: str, ingest_mode: str
) -> Measurement:
    ...


@overload
def _create_or_open_coll(
    cls: Type[Collection[_TDBO]], uri: str, ingest_mode: str
) -> Collection[_TDBO]:
    ...


@typeguard_ignore
def _create_or_open_coll(cls: Type[Any], uri: str, ingest_mode: str) -> Any:
    try:
        thing = cls.open(uri, "w")
    except DoesNotExistError:
        # This is always OK. Make a new one.
        return cls.create(uri)
    # It already exists. Are we resuming?
    if ingest_mode == "resume":
        return thing
    raise SOMAError(f"{uri} already exists")


def _write_dataframe(
    df_uri: str,
    df: pd.DataFrame,
    id_column_name: Optional[str],
    platform_config: Optional[PlatformConfig] = None,
    ingest_mode: IngestMode = "write",
) -> DataFrame:
    s = util.get_start_stamp()
    logging.log_io(None, f"START  WRITING {df_uri}")

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

    try:
        soma_df = factory.open(df_uri, "w", soma_type=DataFrame)
    except DoesNotExistError:
        soma_df = DataFrame.create(
            df_uri, schema=arrow_table.schema, platform_config=platform_config
        )
    else:
        if ingest_mode == "resume":
            storage_ned = _read_nonempty_domain(soma_df)
            dim_range = ((int(df.index.min()), int(df.index.max())),)
            if _chunk_is_contained_in(dim_range, storage_ned):
                logging.log_io(
                    f"Skipped {soma_df.uri}",
                    util.format_elapsed(s, f"SKIPPED {soma_df.uri}"),
                )
                return soma_df
        else:
            raise SOMAError(f"{soma_df.uri} already exists")

    if ingest_mode == "schema_only":
        logging.log_io(
            f"Wrote schema {soma_df.uri}",
            util.format_elapsed(s, f"FINISH WRITING SCHEMA {soma_df.uri}"),
        )
        return soma_df

    soma_df.write(arrow_table)
    logging.log_io(
        f"Wrote   {soma_df.uri}",
        util.format_elapsed(s, f"FINISH WRITING {soma_df.uri}"),
    )
    return soma_df


@typeguard_ignore
def create_from_matrix(
    cls: Type[_NDArr],
    uri: str,
    matrix: Union[Matrix, h5py.Dataset],
    platform_config: Optional[PlatformConfig] = None,
    ingest_mode: IngestMode = "write",
) -> _NDArr:
    """
    Create and populate the ``soma_matrix`` from the contents of ``matrix``.
    """
    # SparseDataset has no ndim but it has a shape
    if len(matrix.shape) != 2:
        raise ValueError(f"expected matrix.shape == 2; got {matrix.shape}")

    s = util.get_start_stamp()
    logging.log_io(None, f"START  WRITING {uri}")

    try:
        soma_ndarray = cls.open(uri, "w", platform_config=platform_config)
    except DoesNotExistError:
        soma_ndarray = cls.create(
            uri,
            type=pa.from_numpy_dtype(matrix.dtype),
            shape=matrix.shape,
            platform_config=platform_config,
        )
    else:
        if ingest_mode != "resume":
            raise SOMAError(f"{soma_ndarray.uri} already exists")

    if ingest_mode == "schema_only":
        logging.log_io(
            f"Wrote schema {soma_ndarray.uri}",
            util.format_elapsed(s, f"FINISH WRITING SCHEMA {soma_ndarray.uri}"),
        )
        return soma_ndarray

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
    elif isinstance(soma_ndarray, SparseNDArray):  # SOMASparseNDArray
        _write_matrix_to_sparseNDArray(
            soma_ndarray,
            matrix,
            tiledb_create_options=TileDBCreateOptions.from_platform_config(
                platform_config
            ),
            ingest_mode=ingest_mode,
        )
    else:
        raise TypeError(f"unknown array type {type(soma_ndarray)}")

    logging.log_io(
        f"Wrote   {soma_ndarray.uri}",
        util.format_elapsed(s, f"FINISH WRITING {soma_ndarray.uri}"),
    )
    return soma_ndarray


def add_X_layer(
    exp: Experiment,
    measurement_name: str,
    X_layer_name: str,
    # E.g. a scipy.csr_matrix from scanpy analysis:
    X_layer_data: Union[Matrix, h5py.Dataset],
    ingest_mode: IngestMode = "write",
    use_relative_uri: Optional[bool] = None,
) -> None:
    """
    This is useful for adding X data, for example from scanpy.pp.normalize_total, scanpy.pp.log1p, etc.

    Use `ingest_mode="resume"` to not error out if the schema already exists.

    [lifecycle: experimental]
    """
    add_matrix_to_collection(
        exp,
        measurement_name,
        "X",
        X_layer_name,
        X_layer_data,
        use_relative_uri=use_relative_uri,
    )


def add_matrix_to_collection(
    exp: Experiment,
    measurement_name: str,
    collection_name: str,
    matrix_name: str,
    # E.g. a scipy.csr_matrix from scanpy analysis:
    matrix_data: Union[Matrix, h5py.Dataset],
    ingest_mode: IngestMode = "write",
    use_relative_uri: Optional[bool] = None,
) -> None:
    """
    This is useful for adding X/obsp/varm/etc data, for example from scanpy.pp.normalize_total,
    scanpy.pp.log1p, etc.

    Use `ingest_mode="resume"` to not error out if the schema already exists.
    """
    with exp.ms[measurement_name] as meas:
        if collection_name in meas:
            coll = cast(Collection[RawHandle], meas[collection_name])
        else:
            coll = _create_or_open_coll(
                Collection, f"{meas.uri}/{collection_name}", ingest_mode
            )
            meas.set(collection_name, coll, use_relative_uri=use_relative_uri)
        with coll:
            uri = f"{coll.uri}/{matrix_name}"
            with create_from_matrix(
                SparseNDArray, uri, matrix_data, ingest_mode=ingest_mode
            ) as sparse_nd_array:
                coll.set(
                    matrix_name, sparse_nd_array, use_relative_uri=use_relative_uri
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
    # TODO: make sure we're not using an old timestamp for this
    storage_ned = None
    if ingest_mode == "resume":
        # This lets us check for already-ingested chunks, when in resume-ingest mode.
        storage_ned = _read_nonempty_domain(soma_ndarray)
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
            chunk_bounds[0] = (
                int(i),
                int(i2 - 1),
            )  # Cast for lint in case np.int64
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


def _read_nonempty_domain(arr: TileDBArray) -> Any:
    try:
        return arr._handle.reader.nonempty_domain()
    except SOMAError:
        # This means that we're open in write-only mode.
        # Reopen the array in read mode.
        pass

    cls = type(arr)
    with cls.open(arr.uri, "r", platform_config=None, context=arr.context) as readarr:
        return readarr._handle.reader.nonempty_domain()


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
    # TODO: make sure we're not using an old timestamp for this
    storage_ned = None
    if ingest_mode == "resume":
        # This lets us check for already-ingested chunks, when in resume-ingest mode.
        # THIS IS A HACK AND ONLY WORKS BECAUSE WE ARE DOING THIS BEFORE ALL WRITES.
        storage_ned = _read_nonempty_domain(soma_ndarray)
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

    [lifecycle: experimental]
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

    [lifecycle: experimental]
    """

    s = util.get_start_stamp()
    logging.log_io(None, "START  Experiment.to_anndata")

    measurement = experiment.ms[measurement_name]

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
        X_csr = conversions.csr_from_tiledb_df(X_mat, nobs, nvar)
        X_dtype = X_csr.dtype
    else:
        raise TypeError(f"Unexpected NDArray type {type(X_data)}")

    obsm = {}
    if "obsm" in measurement:
        for key in measurement.obsm.keys():
            shape = measurement.obsm[key].shape
            if len(shape) != 2:
                raise ValueError(f"expected shape == 2; got {shape}")
            matrix = measurement.obsm[key].read((slice(None),) * len(shape)).to_numpy()
            # The spelling `sp.csr_array` is more idiomatic but doesn't exist until Python 3.8
            obsm[key] = sp.csr_matrix(matrix)

    varm = {}
    if "varm" in measurement:
        for key in measurement.varm.keys():
            shape = measurement.varm[key].shape
            if len(shape) != 2:
                raise ValueError(f"expected shape == 2; got {shape}")
            matrix = measurement.varm[key].read((slice(None),) * len(shape)).to_numpy()
            # The spelling `sp.csr_array` is more idiomatic but doesn't exist until Python 3.8
            varm[key] = sp.csr_matrix(matrix)

    obsp = {}
    if "obsp" in measurement:
        for key in measurement.obsp.keys():
            matrix = measurement.obsp[key].read().tables().concat().to_pandas()
            obsp[key] = conversions.csr_from_tiledb_df(matrix, nobs, nobs)

    varp = {}
    if "varp" in measurement:
        for key in measurement.varp.keys():
            matrix = measurement.varp[key].read().tables().concat().to_pandas()
            varp[key] = conversions.csr_from_tiledb_df(matrix, nvar, nvar)

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
