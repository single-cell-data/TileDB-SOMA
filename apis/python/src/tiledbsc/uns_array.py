from typing import Optional

import numpy as np
import pandas as pd
import scipy.sparse
import tiledb

import tiledbsc.util as util

from .logging import logger
from .tiledb_array import TileDBArray
from .tiledb_group import TileDBGroup


class UnsArray(TileDBArray):
    """
    Holds TileDB storage for an array obtained from the nested `anndata.uns` field.
    """

    # ----------------------------------------------------------------
    def __init__(
        self,
        uri: str,
        name: str,
        parent: Optional[TileDBGroup] = None,
    ):
        """
        See the TileDBObject constructor.
        """
        super().__init__(uri=uri, name=name, parent=parent)

    # ----------------------------------------------------------------
    def from_pandas_dataframe(self, df: pd.DataFrame) -> None:
        """
        Ingests an `UnsArray` into TileDB storage, given a pandas.DataFrame.
        """

        s = util.get_start_stamp()
        logger.info(f"{self._indent}START  WRITING PANDAS.DATAFRAME {self.uri}")

        tiledb.from_pandas(
            uri=self.uri,
            dataframe=df,
            sparse=True,
            allows_duplicates=self._soma_options.allows_duplicates,
            ctx=self._ctx,
        )

        logger.info(
            util.format_elapsed(
                s, f"{self._indent}FINISH WRITING PANDAS.DATAFRAME {self.uri}"
            )
        )

    # ----------------------------------------------------------------
    def _maybe_from_numpyable_object(self, obj) -> bool:
        """
        Nominally for ingest of `uns` nested data from anndata objects. Handles scalar or array values
        -- the former, by wrapping in a 1D array. Maps to TileDB / tiledb.from_numpy storage semantics,
        including UTF-8 handling. Supports dtypes like
        """

        if isinstance(obj, np.ndarray):
            arr = util._to_tiledb_supported_array_type(obj)
            self.from_numpy_ndarray(arr)
            return True

        elif isinstance(obj, list):
            arr = np.asarray(obj)
            self.from_numpy_ndarray(arr)
            return True

        elif "numpy" in str(type(obj)):
            arr = np.asarray([obj])
            arr = util._to_tiledb_supported_array_type(arr)
            self.from_numpy_ndarray(arr)
            return True

        else:
            return False

    # ----------------------------------------------------------------
    def from_numpy_ndarray(self, arr: np.ndarray) -> None:
        """
        Writes a numpy.ndarray to a TileDB array, nominally for ingest of `uns` nested data from anndata
        objects. Mostly tiledb.from_numpy, but with some necessary handling for data with UTF-8 values.
        """

        s = util.get_start_stamp()
        logger.info(f"{self._indent}START  WRITING FROM NUMPY.NDARRAY {self.uri}")

        if "numpy" in str(type(arr)) and str(arr.dtype).startswith("<U"):
            # Note arr.astype('str') does not lead to a successfuly tiledb.from_numpy.
            arr = np.array(arr, dtype="O")

        # overwrite = False
        # if self.exists:
        #     overwrite = True
        #     logger.info(f"{self._indent}Re-using existing array {self.uri}")
        # tiledb.from_numpy(uri=self.uri, array=arr, ctx=self._ctx, overwrite=overwrite)
        # TODO: find the right syntax for update-in-place (tiledb.from_pandas uses `mode`)
        tiledb.from_numpy(uri=self.uri, array=arr, ctx=self._ctx)

        logger.info(
            util.format_elapsed(
                s, f"{self._indent}FINISH WRITING FROM NUMPY.NDARRAY {self.uri}"
            )
        )

    # ----------------------------------------------------------------
    def from_scipy_csr(self, csr: scipy.sparse.csr_matrix) -> None:
        """
        Convert ndarray/(csr|csc)matrix to coo_matrix and ingest into TileDB.

        :param csr: Matrix-like object coercible to a scipy coo_matrix.
        """

        s = util.get_start_stamp()
        logger.info(f"{self._indent}START  WRITING FROM SCIPY.SPARSE.CSR {self.uri}")

        nrows, ncols = csr.shape
        if self.exists():
            logger.info(f"{self._indent}Re-using existing array {self.uri}")
        else:
            self.create_empty_array_for_csr("data", csr.dtype, nrows, ncols)

        self.ingest_data_from_csr(csr)

        logger.info(
            util.format_elapsed(
                s, f"{self._indent}FINISH WRITING FROM SCIPY.SPARSE.CSR {self.uri}"
            )
        )

    # ----------------------------------------------------------------
    def create_empty_array_for_csr(
        self, attr_name: str, matrix_dtype: np.dtype, nrows: int, ncols: int
    ) -> None:
        """
        Create a TileDB 2D sparse array with int dimensions and a single attribute.
        Nominally used for uns data.

        :param matrix_dtype: datatype of the matrix
        :param nrows: number of rows in the matrix
        :param ncols: number of columns in the matrix
        """
        assert isinstance(attr_name, str)

        dom = tiledb.Domain(
            tiledb.Dim(
                name="dim0",
                domain=(0, nrows - 1),
                dtype="int32",
                filters=[tiledb.RleFilter()],
            ),
            tiledb.Dim(
                name="dim1",
                domain=(0, ncols - 1),
                dtype="int32",
                filters=[tiledb.ZstdFilter()],
            ),
            ctx=self._ctx,
        )

        att = tiledb.Attr(
            attr_name, dtype=matrix_dtype, filters=[tiledb.ZstdFilter()], ctx=self._ctx
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
            capacity=100000,
            cell_order="row-major",
            tile_order="col-major",
            ctx=self._ctx,
        )
        tiledb.Array.create(self.uri, sch, ctx=self._ctx)

    # ----------------------------------------------------------------
    def ingest_data_from_csr(self, csr: scipy.sparse.csr_matrix) -> None:
        """
        Convert ndarray/(csr|csc)matrix to coo_matrix and ingest into TileDB.

        :param csr: Matrix-like object coercible to a scipy coo_matrix.
        """

        mat_coo = scipy.sparse.coo_matrix(csr)
        d0 = mat_coo.row
        d1 = mat_coo.col

        with tiledb.open(self.uri, mode="w", ctx=self._ctx) as A:
            A[d0, d1] = mat_coo.data

    # ----------------------------------------------------------------
    # TODO: regardless of which matrix type (numpy.ndarray, scipy.sparse.csr_matrix, etc) was
    # written in, this returns always the same type on readback. Perhaps at write time we can save a
    # metadata tag with the provenance-type of the array, and on readback, try to return the same
    # type.
    def to_matrix(self):
        """
        Reads an uns array from TileDB storage and returns a matrix -- currently, always as numpy.ndarray.
        """
        with tiledb.open(self.uri, ctx=self._ctx) as A:
            df = pd.DataFrame(A[:])
            return df.to_numpy()
