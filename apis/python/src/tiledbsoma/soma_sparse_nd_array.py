import math
import time
from typing import Any, Iterator, List, Literal, Optional, Union

import numpy as np
import pandas as pd
import pyarrow as pa
import scipy.sparse as sp
import tiledb

from . import eta, logging, util, util_arrow, util_scipy
from .soma_collection import SOMACollectionBase
from .tiledb_array import TileDBArray
from .types import Matrix, NTuple, SOMASparseNdCoordinates


class SOMASparseNdArray(TileDBArray):
    """
    Represents ``X`` and others.
    """

    _shape: Optional[NTuple] = None

    def __init__(
        self,
        uri: str,
        *,
        parent: Optional[SOMACollectionBase[Any]] = None,
        ctx: Optional[tiledb.Ctx] = None,
    ):
        """
        Also see the ``TileDBObject`` constructor.
        """

        super().__init__(uri=uri, parent=parent, ctx=ctx)

    @property
    def type(self) -> Literal["SOMASparseNdArray"]:
        return "SOMASparseNdArray"

    def create(
        self,
        type: pa.DataType,
        shape: Union[NTuple, List[int]],
    ) -> "SOMASparseNdArray":
        """
        Create a ``SOMASparseNdArray`` named with the URI.

        :param type: an Arrow type defining the type of each element in the array. If the type is unsupported, an error will be raised.

        :param shape: the length of each domain as a list, e.g., [100, 10]. All lengths must be in the uint64 range.
        """

        # check on shape
        if len(shape) == 0 or any(e <= 0 for e in shape):
            raise ValueError(
                "DenseNdArray shape must be non-zero length tuple of ints > 0"
            )

        if not pa.types.is_primitive(type):
            raise TypeError(
                "Unsupported type - SOMADenseNdArray only supports primtive Arrow types"
            )

        level = self._tiledb_platform_config.string_dim_zstd_level

        dims = []
        for e in shape:
            dim = tiledb.Dim(
                # Use tiledb default names like ``__dim_0``
                domain=(0, e - 1),
                tile=min(e, 2048),  # TODO: PARAMETERIZE,
                dtype=np.uint64,
                filters=[tiledb.ZstdFilter(level=level)],
            )
            dims.append(dim)
        dom = tiledb.Domain(dims, ctx=self._ctx)

        attrs = [
            tiledb.Attr(
                name="data",
                dtype=util_arrow.tiledb_type_from_arrow_type(type),
                filters=[tiledb.ZstdFilter()],
                ctx=self._ctx,
            )
        ]

        # TODO: code-dedupe w/ regard to SOMADenseNdArray. The two creates are
        # almost identical & could share a common parent-class _create() method.
        sch = tiledb.ArraySchema(
            domain=dom,
            attrs=attrs,
            sparse=True,
            allows_duplicates=self._tiledb_platform_config.allows_duplicates,
            offsets_filters=[
                tiledb.DoubleDeltaFilter(),
                tiledb.BitWidthReductionFilter(),
                tiledb.ZstdFilter(),
            ],
            capacity=100000,
            cell_order="row-major",
            tile_order="row-major",
            ctx=self._ctx,
        )

        tiledb.Array.create(self._uri, sch, ctx=self._ctx)

        self._common_create()  # object-type metadata etc

        return self

    @property
    def shape(self) -> NTuple:
        """
        Return length of each dimension, always a list of length ``ndims``
        """
        if self._shape is None:
            with self._tiledb_open() as A:
                self._shape = A.shape
        return self._shape

    def reshape(self, shape: NTuple) -> None:
        raise NotImplementedError("reshape operation not implemented.")

    @property
    def ndims(self) -> int:
        """
        Return number of index columns
        """
        return len(self.shape)

    @property
    def is_sparse(self) -> Literal[True]:
        """
        Returns ``True``.
        """
        return True

    @property
    def nnz(self) -> int:
        """
        Return the number of stored values in the array, including explicitly stored zeros.
        """
        raise NotImplementedError("SOMASparseNdArray.nnz is not implemented.")

    def read_sparse_tensor(
        self,
        coords: SOMASparseNdCoordinates,
        *,
        format: Literal["coo", "csr", "csc"] = "coo",
    ) -> Iterator[Union[pa.SparseCOOTensor, pa.SparseCSCMatrix, pa.SparseCSRMatrix]]:
        """
        Read a use-defined slice of the SparseNdArray and return as an Arrow sparse tensor.

        Parameters
        ----------
        coords : Tuple[Union[int, slice, Tuple[int, ...], List[int], pa.IntegerArray], ...]
            Per-dimension tuple of scalar, slice, sequence of scalar or Arrow IntegerArray
            Arrow arrays currently uninimplemented.

        format - Literal["coo", "csr", "csc"]
            Requested return format:
            * ``coo`` - return an Arrow SparseCOOTensor (default)
            * ``csr`` - return an Arrow SparseCSRMatrix
            * ``csc`` - return an Arrow SparseCSCMatrix


        Returns
        -------
        The requested data in the sparse tensor format specified.
        """

        if format != "coo" and self.ndims != 2:
            raise ValueError(f"Format {format} only supported for 2D SparseNdArray")
        if format not in ("coo", "csr", "csc"):
            raise NotImplementedError("format not implemented")

        with self._tiledb_open("r") as A:
            query = A.query(
                return_arrow=True,
                return_incomplete=True,
            )
            for arrow_tbl in query.df[coords]:
                """
                In PyArrow 9.0.0, there is a bug preventing the creation of "empty"
                (zero element) SparseCOOTensor objects.

                See https://issues.apache.org/jira/browse/ARROW-17933

                Just stop the iteration when we run out of results. The caller must be
                prepared to have a StopIteration, rather than an empty tensor, as the result of
                an empty query.
                """
                if arrow_tbl.num_rows == 0:
                    return

                if format == "coo":

                    coo_data = arrow_tbl.column("data").to_numpy()
                    coo_coords = np.array(
                        [
                            arrow_tbl.column(f"__dim_{n}").to_numpy()
                            for n in range(self.ndims)
                        ]
                    ).T
                    yield pa.SparseCOOTensor.from_numpy(
                        coo_data, coo_coords, shape=A.shape
                    )

                elif format in ("csr", "csc"):
                    # Temporary: as these must be 2D, convert to scipy COO and use
                    # scipy to perform conversions.  C++ reader will be nicer!
                    data = arrow_tbl.column("data").to_numpy()
                    row = arrow_tbl.column("__dim_0").to_numpy()
                    col = arrow_tbl.column("__dim_1").to_numpy()
                    scipy_coo = sp.coo_array((data, (row, col)), shape=A.shape)
                    if format == "csr":
                        yield pa.SparseCSRMatrix.from_scipy(scipy_coo.tocsr())
                    if format == "csc":
                        yield pa.SparseCSCMatrix.from_scipy(scipy_coo.tocsc())

    def read_table(self, coords: SOMASparseNdCoordinates) -> Iterator[pa.Table]:
        """
        Read a user-defined slice of the sparse array and return in COO format
        as an Arrow Table
        """
        with self._tiledb_open("r") as A:
            query = A.query(
                return_arrow=True,
                return_incomplete=True,
            )
            for arrow_tbl in query.df[coords]:
                yield arrow_tbl

    def read_as_pandas(self, coords: SOMASparseNdCoordinates) -> Iterator[pd.DataFrame]:
        """
        Read a user-defined slice of the sparse array and return as a Pandas DataFrame
        containing COO data.
        """
        for arrow_tbl in self.read_table(coords):
            yield arrow_tbl.to_pandas()

    def read_as_pandas_all(
        self, coords: Optional[SOMASparseNdCoordinates] = None
    ) -> pd.DataFrame:
        """
        Return the sparse array as a single Pandas DataFrame containing COO data.
        """
        if coords is None:
            coords = (slice(None),) * self.ndims
        return pd.concat(self.read_as_pandas(coords))

    def write_sparse_tensor(
        self,
        tensor: Union[
            pa.SparseCOOTensor,
            pa.SparseCSRMatrix,
            pa.SparseCSCMatrix,
        ],
    ) -> None:
        """
        Write an Arrow sparse tensor to the SOMASparseNdArray. The coordinates in the Arrow
        SparseTensor will be interpreted as the coordinates to write to.

        Currently supports the _experimental_ Arrow SparseCOOTensor, SparseCSRMatrix and
        SparseCSCMatrix. There is currently no support for Arrow SparseCSFTensor or dense
        Tensor.
        """
        if isinstance(tensor, pa.SparseCOOTensor):
            data, coords = tensor.to_numpy()
            with self._tiledb_open("w") as A:
                A[tuple(c for c in coords.T)] = data
            return

        if isinstance(tensor, (pa.SparseCSCMatrix, pa.SparseCSRMatrix)):
            if self.ndims != 2:
                raise ValueError(
                    f"Unable to write 2D Arrow sparse matrix to {self.ndims}D SOMASparseNdArray"
                )
            sp = tensor.to_scipy().tocoo()
            with self._tiledb_open("w") as A:
                A[sp.row, sp.col] = sp.data
            return

        raise TypeError("Unsuppoted tensor type")

    def write_table(self, arrow_table: pa.Table) -> None:
        """
        Write a COO table, with columns named ``__dim_0``, ..., ``__dim_N`` and ``data``
        to the dense nD array.
        """
        data = arrow_table.column("data").to_numpy()
        coord_tbl = arrow_table.drop(["data"])
        coords = tuple(
            coord_tbl.column(f"__dim_{n}").to_numpy()
            for n in range(coord_tbl.num_columns)
        )
        with self._tiledb_open("w") as A:
            A[coords] = data

    # ---------------
    # TODO: most of this code should be moved into a helper package (eg, soma.io), as
    # is it not part of the core API, and could be built on top of the core API.
    # ---------------
    def from_matrix(self, matrix: Matrix) -> None:
        """
        Imports a matrix -- nominally ``scipy.sparse.csr_matrix`` or ``numpy.ndarray`` -- into a TileDB array which is used for ``X``, ``obsm``, and ``varm`` matrices
        """

        # TODO: Note that scipy.sparse.*_matrix is being deprecated in favor of sparse.*_array.
        # https://docs.scipy.org/doc/scipy/reference/sparse.html

        s = util.get_start_stamp()
        logging.log_io(None, f"{self._indent}START  WRITING")

        if self.exists():
            logging.log_io(None, f"{self._indent}Re-using existing array")
        else:
            self._create_empty_array(
                matrix_dtype=matrix.dtype,
                num_rows=matrix.shape[0],
                num_cols=matrix.shape[1],
            )

        if not self._tiledb_platform_config.write_X_chunked:
            self._ingest_data_whole(matrix)
        elif isinstance(matrix, sp.csr_matrix):
            self._ingest_data_rows_chunked(matrix)
        elif isinstance(matrix, sp.csc_matrix):
            self._ingest_data_cols_chunked(matrix)
        else:
            self._ingest_data_dense_rows_chunked(matrix)

        self._common_create()  # object-type metadata etc

        logging.log_io(
            f"Wrote {self.uri}",
            util.format_elapsed(s, f"{self._indent}FINISH WRITING"),
        )

    # ----------------------------------------------------------------
    def _create_empty_array(
        self, *, matrix_dtype: np.dtype, num_rows: int, num_cols: int
    ) -> None:
        """
        Create a TileDB 2D sparse array with int dimensions and a single attribute.
        """

        dom = tiledb.Domain(
            tiledb.Dim(
                domain=(0, num_rows - 1),
                dtype=np.uint64,
                # TODO: filters=[tiledb.RleFilter()],
            ),
            tiledb.Dim(
                domain=(0, num_cols - 1),
                dtype=np.uint64,
                # TODO: filters=[tiledb.ZstdFilter(level=level)],
            ),
            ctx=self._ctx,
        )

        attrs = tiledb.Attr(
            dtype=matrix_dtype,
            filters=[tiledb.ZstdFilter()],
            ctx=self._ctx,
        )

        sch = tiledb.ArraySchema(
            domain=dom,
            attrs=(attrs,),
            sparse=True,
            allows_duplicates=self._tiledb_platform_config.allows_duplicates,
            offsets_filters=[
                tiledb.DoubleDeltaFilter(),
                tiledb.BitWidthReductionFilter(),
                tiledb.ZstdFilter(),
            ],
            capacity=self._tiledb_platform_config.X_capacity,
            cell_order=self._tiledb_platform_config.X_cell_order,
            tile_order=self._tiledb_platform_config.X_tile_order,
            ctx=self._ctx,
        )

        tiledb.Array.create(self.uri, sch, ctx=self._ctx)

    # ----------------------------------------------------------------
    def _ingest_data_whole(
        self,
        matrix: Matrix,
    ) -> None:
        """
        Convert ``numpy.ndarray``, ``scipy.sparse.csr_matrix``, or ``scipy.sparse.csc_matrix`` to COO matrix and ingest into TileDB.

        :param matrix: Matrix-like object coercible to a scipy COO matrix.
        """

        # TODO: support N-D, not just 2-D
        assert len(matrix.shape) == 2

        mat_coo = sp.coo_matrix(matrix)

        with tiledb.open(self.uri, mode="w", ctx=self._ctx) as A:
            A[mat_coo.row, mat_coo.col] = mat_coo.data

    # ----------------------------------------------------------------
    def _ingest_data_rows_chunked(self, matrix: sp.csr_matrix) -> None:
        """
        Convert csr_matrix to coo_matrix chunkwise and ingest into TileDB.

        :param uri: TileDB URI of the array to be written.
        :param matrix: csr_matrix.
        """

        s = util.get_start_stamp()
        logging.log_io(
            None,
            f"{self._indent}START  ingest",
        )

        nrow = matrix.shape[0]

        eta_tracker = eta.Tracker()
        with tiledb.open(self.uri, mode="w", ctx=self._ctx) as A:

            i = 0
            while i < nrow:
                t1 = time.time()
                # Find a number of CSR rows which will result in a desired nnz for the chunk.
                chunk_size = util_scipy.find_csr_chunk_size(
                    matrix, i, self._tiledb_platform_config.goal_chunk_nnz
                )
                i2 = i + chunk_size

                # Convert the chunk to a COO matrix.
                chunk_coo = matrix[i:i2].tocoo()

                # Python ranges are (lo, hi) with lo inclusive and hi exclusive. But saying that
                # makes us look buggy if we say we're ingesting chunk 0:18 and then 18:32.
                # Instead, print doubly-inclusive lo..hi like 0..17 and 18..31.
                chunk_percent = min(100, 100 * (i2 - 1) / nrow)
                logging.log_io(
                    None,
                    "%sSTART  chunk rows %d..%d of %d (%.3f%%), nnz=%d"
                    % (
                        self._indent,
                        i,
                        i2 - 1,
                        nrow,
                        chunk_percent,
                        chunk_coo.nnz,
                    ),
                )

                # Write a TileDB fragment
                A[chunk_coo.row + i, chunk_coo.col] = chunk_coo.data

                t2 = time.time()
                chunk_seconds = t2 - t1
                eta_seconds = eta_tracker.ingest_and_predict(
                    chunk_percent, chunk_seconds
                )

                if chunk_percent < 100:
                    logging.log_io(
                        "... %7.3f%% done, ETA %s" % (chunk_percent, eta_seconds),
                        "%sFINISH chunk in %.3f seconds, %7.3f%% done, ETA %s"
                        % (self._indent, chunk_seconds, chunk_percent, eta_seconds),
                    )

                i = i2

        logging.log_io(
            None,
            util.format_elapsed(
                s,
                f"{self._indent}FINISH ingest",
            ),
        )

    # This method is very similar to _ingest_data_rows_chunked. The code is largely repeated,
    # and this is intentional. The algorithm here is non-trivial (among the most non-trivial
    # in this package), and adding an abstraction layer would overconfuse it. Here we err
    # on the side of increased readability, at the expense of line-count.
    def _ingest_data_cols_chunked(self, matrix: sp.csc_matrix) -> None:
        """
        Convert csc_matrix to coo_matrix chunkwise and ingest into TileDB.

        :param uri: TileDB URI of the array to be written.
        :param matrix: csc_matrix.
        """

        s = util.get_start_stamp()
        logging.log_io(
            None,
            f"{self._indent}START  ingest",
        )

        ncol = matrix.shape[1]

        eta_tracker = eta.Tracker()
        with tiledb.open(self.uri, mode="w", ctx=self._ctx) as A:

            j = 0
            while j < ncol:
                t1 = time.time()
                # Find a number of CSC columns which will result in a desired nnz for the chunk.
                chunk_size = util_scipy.find_csc_chunk_size(
                    matrix, j, self._tiledb_platform_config.goal_chunk_nnz
                )
                j2 = j + chunk_size

                # Convert the chunk to a COO matrix.
                chunk_coo = matrix[:, j:j2].tocoo()

                # Python ranges are (lo, hi) with lo inclusive and hi exclusive. But saying that
                # makes us look buggy if we say we're ingesting chunk 0:18 and then 18:32.
                # Instead, print doubly-inclusive lo..hi like 0..17 and 18..31.
                chunk_percent = min(100, 100 * (j2 - 1) / ncol)
                logging.log_io(
                    None,
                    "%sSTART  chunk cols %d..%d of %d (%.3f%%), , nnz=%d"
                    % (
                        self._indent,
                        j,
                        j2 - 1,
                        ncol,
                        chunk_percent,
                        chunk_coo.nnz,
                    ),
                )

                # Write a TileDB fragment
                A[chunk_coo.row, chunk_coo.col + j] = chunk_coo.data

                t2 = time.time()
                chunk_seconds = t2 - t1
                eta_seconds = eta_tracker.ingest_and_predict(
                    chunk_percent, chunk_seconds
                )

                if chunk_percent < 100:
                    logging.log_io(
                        "... %7.3f%% done, ETA %s" % (chunk_percent, eta_seconds),
                        "%sFINISH chunk in %.3f seconds, %7.3f%% done, ETA %s"
                        % (self._indent, chunk_seconds, chunk_percent, eta_seconds),
                    )

                j = j2

        logging.log_io(
            None,
            util.format_elapsed(
                s,
                f"{self._indent}FINISH ingest",
            ),
        )

    # This method is very similar to _ingest_data_rows_chunked. The code is largely repeated,
    # and this is intentional. The algorithm here is non-trivial (among the most non-trivial
    # in this package), and adding an abstraction layer would overconfuse it. Here we err
    # on the side of increased readability, at the expense of line-count.
    def _ingest_data_dense_rows_chunked(
        self,
        matrix: np.ndarray,
    ) -> None:
        """
        Convert dense matrix to coo_matrix chunkwise and ingest into TileDB.

        :param uri: TileDB URI of the array to be written.
        :param matrix: dense matrix.
        """

        nrow, ncol = matrix.shape

        s = util.get_start_stamp()
        logging.log_io(
            None,
            f"{self._indent}START  ingest",
        )

        eta_tracker = eta.Tracker()
        with tiledb.open(self.uri, mode="w", ctx=self._ctx) as A:

            i = 0
            while i < nrow:
                t1 = time.time()
                # Find a number of dense rows which will result in a desired nnz for the chunk,
                # rounding up to the nearest integer. Example: goal_chunk_nnz is 120. ncol is 50;
                # 120/50 rounds up to 3; take 3 rows per chunk.
                chunk_size = int(
                    math.ceil(self._tiledb_platform_config.goal_chunk_nnz / ncol)
                )
                i2 = i + chunk_size

                # Convert the chunk to a COO matrix.
                chunk = matrix[i:i2]
                chunk_coo = sp.csr_matrix(chunk).tocoo()

                # Python ranges are (lo, hi) with lo inclusive and hi exclusive. But saying that
                # makes us look buggy if we say we're ingesting chunk 0:18 and then 18:32.
                # Instead, print doubly-inclusive lo..hi like 0..17 and 18..31.
                chunk_percent = min(100, 100 * (i2 - 1) / nrow)
                logging.log_io(
                    None,
                    "%sSTART  chunk rows %d..%d of %d (%.3f%%), nnz=%d"
                    % (
                        self._indent,
                        i,
                        i2 - 1,
                        nrow,
                        chunk_percent,
                        chunk_coo.nnz,
                    ),
                )

                # Write a TileDB fragment
                A[chunk_coo.row + i, chunk_coo.col] = chunk_coo.data

                t2 = time.time()
                chunk_seconds = t2 - t1
                eta_seconds = eta_tracker.ingest_and_predict(
                    chunk_percent, chunk_seconds
                )

                if chunk_percent < 100:
                    logging.log_io(
                        "... %7.3f%% done, ETA %s" % (chunk_percent, eta_seconds),
                        "%sFINISH chunk in %.3f seconds, %7.3f%% done, ETA %s"
                        % (self._indent, chunk_seconds, chunk_percent, eta_seconds),
                    )

                i = i2

        logging.log_io(
            None,
            util.format_elapsed(
                s,
                f"{self._indent}FINISH ingest",
            ),
        )
