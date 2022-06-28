import math
import time
from typing import Optional, Tuple, Union

import numpy as np
import pandas as pd
import scipy.sparse
import tiledb

import tiledbsc.util as util

from .annotation_dataframe import AnnotationDataFrame
from .logging import logger
from .tiledb_array import TileDBArray
from .tiledb_group import TileDBGroup


class AssayMatrix(TileDBArray):
    """
    Wraps a TileDB sparse array with two string dimensions.
    Used for `X`, `raw.X`, `obsp` elements, and `varp` elements.
    """

    row_dim_name: str  # obs_id for X, obs_id_i for obsp; var_id_i for varp
    col_dim_name: str  # var_id for X, obs_id_j for obsp; var_id_j for varp
    attr_name: str
    row_dataframe: AnnotationDataFrame
    col_dataframe: AnnotationDataFrame

    # ----------------------------------------------------------------
    def __init__(
        self,
        uri: str,
        name: str,
        row_dim_name: str,
        col_dim_name: str,
        row_dataframe: AnnotationDataFrame,  # Nominally a reference to soma.obs
        col_dataframe: AnnotationDataFrame,  # Nominally a reference to soma.var
        parent: Optional[TileDBGroup] = None,
    ):
        """
        See the TileDBObject constructor.

        The `row_dataframe` and `col_dataframe` are nominally:

        * `soma.obs` and `soma.var`, for `soma.X["data"]`
        * `soma.obs` and `soma.raw.var`, for `soma.raw.X["data"]`
        * `soma.obs` and `soma.obs`, for `soma.obsp` elements
        * `soma.var` and `soma.var`, for `soma.obsp` elements

        References to these objects are kept solely for obtaining dim labels for metadata
        acquisition at runtime (e.g. shape). We retain references to these objects, rather
        than taking in actualized ID-lists here in the constructor, for two reasons:

        * We need to be able to set up a SOMA to write to, before it's been populated.
        * For reading from an already-populated SOMA, we wish to avoid cache-coherency issues.
        """
        super().__init__(uri=uri, name=name, parent=parent)

        self.row_dim_name = row_dim_name
        self.col_dim_name = col_dim_name
        self.attr_name = "value"
        self.row_dataframe = row_dataframe
        self.col_dataframe = col_dataframe

    # ----------------------------------------------------------------
    def shape(self) -> Tuple[int, int]:
        """
        Returns a tuple with the number of rows and number of columns of the `AssayMatrix`.
        In TileDB storage, these are string-indexed sparse arrays for which no `.shape()` exists,
        but, we draw from the appropriate `obs`, `var`, `raw/var`, etc. as appropriate for a given matrix.

        Note: currently implemented via data scan -- will be optimized for TileDB core 2.10.
        """
        with self._open():
            # These TileDB arrays are string-dimensioned sparse arrays so there is no '.shape'.
            # Instead we compute it ourselves.  See also:
            num_rows = self.row_dataframe.shape()[0]
            num_cols = self.col_dataframe.shape()[0]
            return (num_rows, num_cols)

    # ----------------------------------------------------------------
    def dim_select(self, obs_ids, var_ids) -> pd.DataFrame:
        """
        Selects a slice out of the matrix with specified `obs_ids` and/or `var_ids`.
        Either or both of the ID lists may be `None`, meaning, do not subselect along
        that dimension. If both ID lists are `None`, the entire matrix is returned.
        """
        with tiledb.open(self.uri, ctx=self._ctx) as A:
            if obs_ids is None:
                if var_ids is None:
                    df = A.df[:, :]
                else:
                    df = A.df[:, var_ids]
            else:
                if var_ids is None:
                    df = A.df[obs_ids, :]
                else:
                    df = A.df[obs_ids, var_ids]
        df.set_index([self.row_dim_name, self.col_dim_name], inplace=True)
        return df

    # ----------------------------------------------------------------
    def df(self, obs_ids=None, var_ids=None) -> pd.DataFrame:
        """
        Keystroke-saving alias for `.dim_select()`. If either of `obs_ids` or `var_ids`
        are provided, they're used to subselect; if not, the entire dataframe is returned.
        """
        return self.dim_select(obs_ids, var_ids)

    # ----------------------------------------------------------------
    def csr(self, obs_ids=None, var_ids=None) -> scipy.sparse.csr_matrix:
        """
        Like `.df()` but returns results in `scipy.sparse.csr_matrix` format.
        """
        return self._csr_or_csc("csr", obs_ids, var_ids)

    def csc(self, obs_ids=None, var_ids=None) -> scipy.sparse.csc_matrix:
        """
        Like `.df()` but returns results in `scipy.sparse.csc_matrix` format.
        """
        return self._csr_or_csc("csc", obs_ids, var_ids)

    def _csr_or_csc(
        self, which: str, obs_ids=None, var_ids=None
    ) -> Union[scipy.sparse.csr_matrix, scipy.sparse.csc_matrix]:
        """
        Helper method for `csr` and `csc`.
        """
        assert which in ("csr", "csc")
        df = self.dim_select(obs_ids, var_ids)
        if obs_ids is None:
            obs_ids = self.row_dataframe.ids()
        if var_ids is None:
            var_ids = self.col_dataframe.ids()
        return util.X_and_ids_to_sparse_matrix(
            df,
            self.row_dim_name,
            self.col_dim_name,
            self.attr_name,
            obs_ids,
            var_ids,
            which,
        )

    # ----------------------------------------------------------------
    def from_matrix_and_dim_values(self, matrix, row_names, col_names) -> None:
        """
        Imports a matrix -- nominally `scipy.sparse.csr_matrix` or `numpy.ndarray` -- into a TileDB
        array which is used for `X`, `raw.X`, `obsp` members, and `varp` members.

        The `row_names` and `col_names` are row and column labels for the matrix; the matrix may be
        `scipy.sparse.csr_matrix`, `scipy.sparse.csc_matrix`, `numpy.ndarray`, etc.
        For ingest from `AnnData`, these should be `ann.obs_names` and `ann.var_names`.
        """

        s = util.get_start_stamp()
        logger.info(f"{self._indent}START  WRITING {self.uri}")

        assert len(row_names) == matrix.shape[0]
        assert len(col_names) == matrix.shape[1]

        # Following Pythonic practice, the row_names and col_names can be all manner of things:
        # pandas.core.indexes.base.Index, numpy.ndarray, list of string, etc. However, we do have
        # one requirement: that they be addressable via multi-index like `row_names[[0,1,2]]`.
        if isinstance(row_names, list):
            row_names = np.asarray(row_names)
        if isinstance(col_names, list):
            col_names = np.asarray(col_names)

        if self.exists():
            logger.info(f"{self._indent}Re-using existing array {self.uri}")
        else:
            self._create_empty_array(matrix_dtype=matrix.dtype)

        self._set_object_type_metadata()

        self._ingest_data(matrix, row_names, col_names)
        logger.info(util.format_elapsed(s, f"{self._indent}FINISH WRITING {self.uri}"))

    # ----------------------------------------------------------------
    def _create_empty_array(self, matrix_dtype: np.dtype) -> None:
        """
        Create a TileDB 2D sparse array with string dimensions and a single attribute.
        """

        level = self._soma_options.string_dim_zstd_level
        dom = tiledb.Domain(
            tiledb.Dim(
                name=self.row_dim_name,
                domain=(None, None),
                dtype="ascii",
                filters=[tiledb.RleFilter()],
            ),
            tiledb.Dim(
                name=self.col_dim_name,
                domain=(None, None),
                dtype="ascii",
                filters=[tiledb.ZstdFilter(level=level)],
            ),
            ctx=self._ctx,
        )

        att = tiledb.Attr(
            self.attr_name,
            dtype=matrix_dtype,
            filters=[tiledb.ZstdFilter()],
            ctx=self._ctx,
        )

        sch = tiledb.ArraySchema(
            domain=dom,
            attrs=(att,),
            sparse=True,
            allows_duplicates=self._soma_options.allows_duplicates,
            offsets_filters=[
                tiledb.DoubleDeltaFilter(),
                tiledb.BitWidthReductionFilter(),
                tiledb.ZstdFilter(),
            ],
            capacity=self._soma_options.X_capacity,
            cell_order=self._soma_options.X_cell_order,
            tile_order=self._soma_options.X_tile_order,
            ctx=self._ctx,
        )

        tiledb.Array.create(self.uri, sch, ctx=self._ctx)

    # ----------------------------------------------------------------
    def _ingest_data(self, matrix, row_names, col_names) -> None:
        if self._soma_options.write_X_chunked:
            if isinstance(matrix, scipy.sparse.csr_matrix):
                self.ingest_data_rows_chunked(matrix, row_names, col_names)
            elif isinstance(matrix, scipy.sparse.csc_matrix):
                self.ingest_data_cols_chunked(matrix, row_names, col_names)
            else:
                self.ingest_data_dense_rows_chunked(matrix, row_names, col_names)
        else:
            self.ingest_data_whole(matrix, row_names, col_names)

    # ----------------------------------------------------------------
    def ingest_data_whole(self, matrix, row_names, col_names) -> None:
        """
        Convert `numpy.ndarray`, `scipy.sparse.csr_matrix`, or `scipy.sparse.csc_matrix` to COO matrix and ingest into TileDB.

        :param matrix: Matrix-like object coercible to a scipy COO matrix.
        :param row_names: List of row names.
        :param col_names: List of column names.
        """

        assert len(row_names) == matrix.shape[0]
        assert len(col_names) == matrix.shape[1]

        mat_coo = scipy.sparse.coo_matrix(matrix)
        d0 = row_names[mat_coo.row]
        d1 = col_names[mat_coo.col]

        with tiledb.open(self.uri, mode="w", ctx=self._ctx) as A:
            A[d0, d1] = mat_coo.data

    # ----------------------------------------------------------------
    # Example: suppose this 4x3 is to be written in two chunks of two rows each
    # but written in sorted order.
    #
    # Original     Sorted     Permutation
    #  data       row names
    #
    #   X Y Z
    # C 0 1 2      A            1
    # A 4 0 5      B            2
    # B 7 0 0      C            0
    # D 0 8 9      D            3
    #
    # First chunk:
    # * Row indices 0,1 map to permutation indices 1,2
    # * i,i2 are 0,2
    # * chunk_coo is original matrix rows 1,2
    # * chunk_coo.row is [0,1]
    # * chunk_coo.row + i is [0,1]
    # * sorted_row_names: ['A', 'B']
    #
    # Second chunk:
    # * Row indices 2,3 map to permutation indices 0,3
    # * i,i2 are 2,4
    # * chunk_coo is original matrix rows 0,3
    # * chunk_coo.row is [0,1]
    # * chunk_coo.row + i is [2,3]
    # * sorted_row_names: ['C', 'D']
    #
    # See README-csr-ingest.md for important information of using this ingestor.
    # ----------------------------------------------------------------

    def ingest_data_rows_chunked(self, matrix, row_names, col_names) -> None:
        """
        Convert csr_matrix to coo_matrix chunkwise and ingest into TileDB.

        :param uri: TileDB URI of the array to be written.
        :param matrix: csr_matrix.
        :param row_names: List of row names.
        :param col_names: List of column names.
        """

        assert len(row_names) == matrix.shape[0]
        assert len(col_names) == matrix.shape[1]

        # Sort the row names so we can write chunks indexed by sorted string keys.  This will lead
        # to efficient TileDB fragments in the sparse array indexed by these string keys.
        #
        # Key note: only the _obs labels_ are being sorted, and along with them come permutation
        # indices for accessing the CSR matrix via cursor-indirection -- e.g. csr row 28 is accessed as
        # csr row permuation[28] -- the CSR matrix itself isn't sorted in bulk.
        sorted_row_names, permutation = util._get_sort_and_permutation(list(row_names))
        # Using numpy we can index this with a list of indices, which a plain Python list doesn't support.
        sorted_row_names = np.asarray(sorted_row_names)

        s = util.get_start_stamp()
        logger.info(f"{self._indent}START  __ingest_coo_data_string_dims_rows_chunked")

        eta_tracker = util.ETATracker()
        with tiledb.open(self.uri, mode="w", ctx=self._ctx) as A:
            nrow = len(sorted_row_names)

            i = 0
            while i < nrow:
                t1 = time.time()
                # Find a number of CSR rows which will result in a desired nnz for the chunk.
                chunk_size = util._find_csr_chunk_size(
                    matrix, permutation, i, self._soma_options.goal_chunk_nnz
                )
                i2 = i + chunk_size

                # Convert the chunk to a COO matrix.
                chunk_coo = matrix[permutation[i:i2]].tocoo()

                # Write the chunk-COO to TileDB.
                d0 = sorted_row_names[chunk_coo.row + i]
                d1 = col_names[chunk_coo.col]

                if len(d0) == 0:
                    i = i2
                    continue

                # Python ranges are (lo, hi) with lo inclusive and hi exclusive. But saying that
                # makes us look buggy if we say we're ingesting chunk 0:18 and then 18:32.
                # Instead, print doubly-inclusive lo..hi like 0..17 and 18..31.
                chunk_percent = 100 * (i2 - 1) / nrow
                logger.info(
                    "%sSTART  chunk rows %d..%d of %d (%.3f%%), obs_ids %s..%s, nnz=%d"
                    % (
                        self._indent,
                        i,
                        i2 - 1,
                        nrow,
                        chunk_percent,
                        d0[0],
                        d0[-1],
                        chunk_coo.nnz,
                    )
                )

                # Write a TileDB fragment
                A[d0, d1] = chunk_coo.data

                t2 = time.time()
                chunk_seconds = t2 - t1
                eta = eta_tracker.ingest_and_predict(chunk_percent, chunk_seconds)

                logger.info(
                    "%sFINISH chunk in %.3f seconds, %7.3f%% done, ETA %s"
                    % (self._indent, chunk_seconds, chunk_percent, eta)
                )

                i = i2

        logger.info(
            util.format_elapsed(
                s,
                f"{self._indent}FINISH __ingest_coo_data_string_dims_rows_chunked",
            )
        )

    # This method is very similar to ingest_data_rows_chunked. The code is largely repeated,
    # and this is intentional. The algorithm here is non-trivial (among the most non-trivial
    # in this package), and adding an abstraction layer would overconfuse it. Here we err
    # on the side of increased readability, at the expense of line-count.
    def ingest_data_cols_chunked(self, matrix, row_names, col_names) -> None:
        """
        Convert csc_matrix to coo_matrix chunkwise and ingest into TileDB.

        :param uri: TileDB URI of the array to be written.
        :param matrix: csc_matrix.
        :param row_names: List of row names.
        :param col_names: List of column names.
        """

        assert len(row_names) == matrix.shape[0]
        assert len(col_names) == matrix.shape[1]

        # Sort the column names so we can write chunks indexed by sorted string keys.  This will lead
        # to efficient TileDB fragments in the sparse array indexed by these string keys.
        #
        # Key note: only the _var labels_ are being sorted, and along with them come permutation
        # indices for accessing the CSC matrix via cursor-indirection -- e.g. csc column 28 is
        # accessed as csc column permuation[28] -- the CSC matrix itself isn't sorted in bulk.
        sorted_col_names, permutation = util._get_sort_and_permutation(list(col_names))
        # Using numpy we can index this with a list of indices, which a plain Python list doesn't support.
        sorted_col_names = np.asarray(sorted_col_names)

        s = util.get_start_stamp()
        logger.info(f"{self._indent}START  __ingest_coo_data_string_dims_cols_chunked")

        eta_tracker = util.ETATracker()
        with tiledb.open(self.uri, mode="w", ctx=self._ctx) as A:
            ncol = len(sorted_col_names)

            j = 0
            while j < ncol:
                t1 = time.time()
                # Find a number of CSC columns which will result in a desired nnz for the chunk.
                chunk_size = util._find_csc_chunk_size(
                    matrix, permutation, j, self._soma_options.goal_chunk_nnz
                )
                j2 = j + chunk_size

                # Convert the chunk to a COO matrix.
                chunk_coo = matrix[:, permutation[j:j2]].tocoo()

                # Write the chunk-COO to TileDB.
                d0 = row_names[chunk_coo.row]
                d1 = sorted_col_names[chunk_coo.col + j]

                if len(d1) == 0:
                    j = j2
                    continue

                # Python ranges are (lo, hi) with lo inclusive and hi exclusive. But saying that
                # makes us look buggy if we say we're ingesting chunk 0:18 and then 18:32.
                # Instead, print doubly-inclusive lo..hi like 0..17 and 18..31.
                chunk_percent = 100 * (j2 - 1) / ncol
                logger.info(
                    "%sSTART  chunk rows %d..%d of %d (%.3f%%), var_ids %s..%s, nnz=%d"
                    % (
                        self._indent,
                        j,
                        j2 - 1,
                        ncol,
                        chunk_percent,
                        d1[0],
                        d1[-1],
                        chunk_coo.nnz,
                    )
                )

                # Write a TileDB fragment
                A[d0, d1] = chunk_coo.data

                t2 = time.time()
                chunk_seconds = t2 - t1
                eta = eta_tracker.ingest_and_predict(chunk_percent, chunk_seconds)

                logger.info(
                    "%sFINISH chunk in %.3f seconds, %7.3f%% done, ETA %s"
                    % (self._indent, chunk_seconds, chunk_percent, eta)
                )

                j = j2

        logger.info(
            util.format_elapsed(
                s,
                f"{self._indent}FINISH __ingest_coo_data_string_dims_rows_chunked",
            )
        )

    # This method is very similar to ingest_data_rows_chunked. The code is largely repeated,
    # and this is intentional. The algorithm here is non-trivial (among the most non-trivial
    # in this package), and adding an abstraction layer would overconfuse it. Here we err
    # on the side of increased readability, at the expense of line-count.
    def ingest_data_dense_rows_chunked(self, matrix, row_names, col_names) -> None:
        """
        Convert dense matrix to coo_matrix chunkwise and ingest into TileDB.

        :param uri: TileDB URI of the array to be written.
        :param matrix: dense matrix.
        :param row_names: List of row names.
        :param col_names: List of column names.
        """

        assert len(row_names) == matrix.shape[0]
        assert len(col_names) == matrix.shape[1]

        # Sort the row names so we can write chunks indexed by sorted string keys.  This will lead
        # to efficient TileDB fragments in the sparse array indexed by these string keys.
        #
        # Key note: only the _obs labels_ are being sorted, and along with them come permutation
        # indices for accessing the dense matrix via cursor-indirection -- e.g. dense row 28 is accessed as
        # dense row permuation[28] -- the dense matrix itself isn't sorted in bulk.
        sorted_row_names, permutation = util._get_sort_and_permutation(list(row_names))
        # Using numpy we can index this with a list of indices, which a plain Python list doesn't support.
        sorted_row_names = np.asarray(sorted_row_names)

        s = util.get_start_stamp()
        logger.info(
            f"{self._indent}START  __ingest_coo_data_string_dims_dense_rows_chunked"
        )

        eta_tracker = util.ETATracker()
        with tiledb.open(self.uri, mode="w", ctx=self._ctx) as A:
            nrow = len(sorted_row_names)
            ncol = len(col_names)

            i = 0
            while i < nrow:
                t1 = time.time()
                # Find a number of dense rows which will result in a desired nnz for the chunk,
                # rounding up to the nearest integer. Example: goal_chunk_nnz is 120. ncol is 50;
                # 120/50 rounds up to 3; take 3 rows per chunk.
                chunk_size = int(math.ceil(self._soma_options.goal_chunk_nnz / ncol))
                i2 = i + chunk_size

                # Convert the chunk to a COO matrix.
                chunk = matrix[permutation[i:i2]]
                chunk_coo = scipy.sparse.csr_matrix(chunk).tocoo()

                # Write the chunk-COO to TileDB.
                d0 = sorted_row_names[chunk_coo.row + i]
                d1 = col_names[chunk_coo.col]

                if len(d0) == 0:
                    i = i2
                    continue

                # Python ranges are (lo, hi) with lo inclusive and hi exclusive. But saying that
                # makes us look buggy if we say we're ingesting chunk 0:18 and then 18:32.
                # Instead, print doubly-inclusive lo..hi like 0..17 and 18..31.
                chunk_percent = 100 * (i2 - 1) / nrow
                logger.info(
                    "%sSTART  chunk rows %d..%d of %d (%.3f%%), obs_ids %s..%s, nnz=%d"
                    % (
                        self._indent,
                        i,
                        i2 - 1,
                        nrow,
                        chunk_percent,
                        d0[0],
                        d0[-1],
                        chunk_coo.nnz,
                    )
                )

                # Write a TileDB fragment
                A[d0, d1] = chunk_coo.data

                t2 = time.time()
                chunk_seconds = t2 - t1
                eta = eta_tracker.ingest_and_predict(chunk_percent, chunk_seconds)

                logger.info(
                    "%sFINISH chunk in %.3f seconds, %7.3f%% done, ETA %s"
                    % (self._indent, chunk_seconds, chunk_percent, eta)
                )

                i = i2

        logger.info(
            util.format_elapsed(
                s,
                f"{self._indent}FINISH __ingest_coo_data_string_dims_dense_rows_chunked",
            )
        )

    # ----------------------------------------------------------------
    def to_csr_matrix(self, row_labels, col_labels) -> scipy.sparse.csr_matrix:
        """
        Reads the TileDB array storage for the storage and returns a sparse CSR matrix.  The
        row/columns labels should be `obs,var` labels if the `AssayMatrix` is `X`, or `obs,obs` labels if
        the `AssayMatrix` is `obsp`, or `var,var` labels if the `AssayMatrix` is `varp`.
        Note in all cases that TileDB will have sorted the row and column labels; they won't
        be in the same order as they were in any anndata object which was used to create the
        TileDB storage.
        """

        s = util.get_start_stamp()
        logger.info(f"{self._indent}START  read {self.uri}")

        csr = self.csr()

        logger.info(util.format_elapsed(s, f"{self._indent}FINISH read {self.uri}"))

        return csr
