import math
import time
from typing import List, Optional, Sequence, Union

import numpy as np
import pandas as pd
import pyarrow as pa
import scipy.sparse as sp
import tiledb

import tiledbsc.v1.eta as eta
import tiledbsc.v1.logging as logging
import tiledbsc.v1.util as util
import tiledbsc.v1.util_arrow as util_arrow

from .soma_collection import SOMACollection
from .tiledb_array import TileDBArray
from .types import Matrix, NTuple


class SOMASparseNdArray(TileDBArray):
    """
    Represents ``X`` and others.
    """

    _shape: Optional[NTuple] = None

    def __init__(
        self,
        uri: str,
        *,
        name: Optional[str] = None,
        parent: Optional[SOMACollection] = None,
    ):
        """
        Also see the :class:`TileDBObject` constructor.
        """

        super().__init__(uri=uri, name=name, parent=parent)

    def create(
        self,
        type: pa.DataType,
        shape: Union[NTuple, List[int]],
    ) -> None:
        """
        Create a SOMASparseNdArray named with the URI.

        :param type: an Arrow type defining the type of each element in the array. If the type is
        unsupported, an error will be raised.

        :param shape: the length of each domain as a list, e.g., [100, 10]. All lengths must be in
        the uint64 range.
        """

        # checks on shape
        assert len(shape) > 0
        for e in shape:
            assert e > 0

        level = self._tiledb_platform_config.string_dim_zstd_level

        dims = []
        for e in shape:
            dim = tiledb.Dim(
                # Use tiledb default names like `__dim_0`
                domain=(0, e - 1),
                tile=min(e, 2048),  # TODO: PARAMETERIZE
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

    def __repr__(self) -> str:
        """
        Default display of `SOMASparseNdArray`.
        """
        return "\n".join(self._repr_aux())

    def _repr_aux(self) -> List[str]:
        lines = [
            self.get_name()
            + " "
            + self.__class__.__name__
            + " "
            + str(self.get_shape())
        ]
        return lines

    def get_shape(self) -> NTuple:
        """
        Return length of each dimension, always a list of length ``ndims``
        """
        if self._shape is None:
            with self._tiledb_open() as A:
                self._shape = A.shape
        return self._shape

    def get_ndims(self) -> int:
        """
        Return number of index columns
        """
        return len(self.get_shape())

    # TODO
    #    def get_schema(self) -> Arrow.Schema:
    #        """
    #        Return data schema, in the form of an Arrow Schema
    #        """

    def get_is_sparse(self) -> bool:
        """
        Returns ``True``.
        """
        return True

    # TODO
    #    def get_nnz(self) -> wint:
    #        """
    #        Return the number of non-zero values in the array
    #        """
    #        return 999

    # TODO
    #    def read():
    #        """
    #        Read a slice of data from the SOMASparseNdArray
    #        """

    # ### Operation: read()
    #
    # Read a user-specified subset of the object, and return as one or more Arrow.SparseTensor.
    #
    # Summary:
    #
    # ```
    # read(
    #     [slice, ...],
    #     partitions,
    #     result_order
    # ) -> delayed iterator over Arrow.SparseTensor
    # ```
    #
    # - slice - per-dimension slice, expressed as a scalar, a range, or a list of both.
    # - partitions - an optional [`SOMAReadPartitions`](#SOMAReadPartitions) hint to indicate how
    #   results should be organized.
    # - result_order - order of read results. Can be one of row-major, column-major and unordered.
    #
    # The `read` operation will return a language-specific iterator over one or more Arrow SparseTensor
    # objects, allowing the incremental processing of results larger than available memory. The actual
    # iterator used is delegated to language-specific SOMA specs.

    def write(
        self,
        # TODO: define a compound type for this in types.py
        tensor: Union[pa.SparseCOOTensor, pa.SparseCSFTensor],
    ) -> None:
        """
        Write an Arrow.Tensor to the persistent object. As duplicate index values are not allowed, index
        values already present in the object are overwritten and new index values are added.

        :param values: an Arrow.SparseTensor containing values to be written. The type of elements in `values`
        must match the type of the `SOMASparseNdArray`.
        """

        # Example caller-side:
        #
        # tensor = pa.SparseCOOTensor.from_numpy(
        #     data=np.asarray([4,5,6]),
        #     coords=[[1,2], [3,4], [5,6]],
        #     shape=(nr, nc),
        # )
        #
        # Here on the callee side: need to figure out how to extrat the coords separate from the
        # data since tiledb write needs `A[coords] = data`.

        # >>> tensor.data
        # AttributeError: 'pyarrow.lib.SparseCOOTensor' object has no attribute 'data'

        # >>> tensor.coords
        # AttributeError: 'pyarrow.lib.SparseCOOTensor' object has no attribute 'coords'

        # :headdesk:

        nt = tensor.to_numpy()
        assert len(nt) == 2
        coords = nt[1]
        data = nt[0]

        # The coords come in as a list of [i,j] pairs. We need a list of i's and a list of j's.
        # TODO: for now, only support 2D matrices. We need to complexify this code to handle n-d arrays.
        assert len(coords) > 0
        assert len(coords[0]) == 2
        icoords = coords[:, 0]
        jcoords = coords[:, 1]

        with self._tiledb_open("w") as A:
            A[icoords, jcoords] = data

    # ================================================================
    # ================================================================
    # ================================================================
    def to_pandas(
        self,
        *,
        row_ids: Optional[Sequence[int]] = None,
        col_ids: Optional[Sequence[int]] = None,
        set_index: Optional[bool] = False,
    ) -> pd.DataFrame:
        """
        TODO: comment
        """
        with self._tiledb_open() as A:
            if row_ids is None:
                if col_ids is None:
                    df = A.df[:, :]
                else:
                    df = A.df[:, col_ids]
            else:
                if col_ids is None:
                    df = A.df[row_ids, :]
                else:
                    df = A.df[row_ids, col_ids]
        # Make this opt-in only.  For large arrays, this df.set_index is time-consuming
        # so we should not do it without direction.
        if set_index:
            # TODO: get these from the schema while self._tiledb_open.
            df.set_index(["__dim_0", "__dim_1"], inplace=True)
        return df

    # ----------------------------------------------------------------
    def from_matrix(self, matrix: Matrix) -> None:
        """
        Imports a matrix -- nominally `scipy.sparse.csr_matrix` or `numpy.ndarray` -- into a TileDB
        array which is used for `X`, `obsm`, and `varm` matrices
        """

        s = util.get_start_stamp()
        logging.log_io(None, f"{self._indent}START  WRITING {self.get_uri()}")

        if self.exists():
            logging.log_io(
                None, f"{self._indent}Re-using existing array {self.get_uri()}"
            )
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
            f"Wrote {self._nested_name}",
            util.format_elapsed(s, f"{self._indent}FINISH WRITING {self.get_uri()}"),
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

        tiledb.Array.create(self.get_uri(), sch, ctx=self._ctx)

    # ----------------------------------------------------------------
    def _ingest_data_whole(
        self,
        matrix: Matrix,
    ) -> None:
        """
        Convert `numpy.ndarray`, `scipy.sparse.csr_matrix`, or `scipy.sparse.csc_matrix`
        to COO matrix and ingest into TileDB.

        :param matrix: Matrix-like object coercible to a scipy COO matrix.
        """

        # TODO: support N-D, not just 2-D
        assert len(matrix.shape) == 2

        mat_coo = sp.coo_matrix(matrix)

        with tiledb.open(self.get_uri(), mode="w", ctx=self._ctx) as A:
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
        with tiledb.open(self.get_uri(), mode="w", ctx=self._ctx) as A:

            i = 0
            while i < nrow:
                t1 = time.time()
                # Find a number of CSR rows which will result in a desired nnz for the chunk.
                chunk_size = util._find_csr_chunk_size(
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
                        "... %s %7.3f%% done, ETA %s"
                        % (self._nested_name, chunk_percent, eta_seconds),
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
        with tiledb.open(self.get_uri(), mode="w", ctx=self._ctx) as A:

            j = 0
            while j < ncol:
                t1 = time.time()
                # Find a number of CSC columns which will result in a desired nnz for the chunk.
                chunk_size = util._find_csc_chunk_size(
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
                        "... %s %7.3f%% done, ETA %s"
                        % (self._nested_name, chunk_percent, eta_seconds),
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
        with tiledb.open(self.get_uri(), mode="w", ctx=self._ctx) as A:

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
                        "... %s %7.3f%% done, ETA %s"
                        % (self._nested_name, chunk_percent, eta_seconds),
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
