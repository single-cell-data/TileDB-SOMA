import collections.abc
from typing import Any, Iterator, List, Literal, Optional, Union, cast

import numpy as np
import pandas as pd
import pyarrow as pa
import scipy.sparse as sp
import tiledb

# This package's pybind11 code
import tiledbsoma.libtiledbsoma as clib

from . import tiledb_platform_config as tdbpc
from . import util, util_arrow
from .collection import CollectionBase
from .tiledb_array import TileDBArray
from .types import NTuple, PlatformConfig, SparseNdCoordinates


class SparseNDArray(TileDBArray):
    """
    Represents ``X`` and others.
    """

    def __init__(
        self,
        uri: str,
        *,
        parent: Optional[CollectionBase[Any]] = None,
        tiledb_platform_config: Optional[tdbpc.TileDBPlatformConfig] = None,
        ctx: Optional[tiledb.Ctx] = None,
    ):
        """
        Also see the ``TileDBObject`` constructor.
        """

        super().__init__(
            uri=uri,
            parent=parent,
            tiledb_platform_config=tiledb_platform_config,
            ctx=ctx,
        )

    @property
    def soma_type(self) -> Literal["SOMASparseNDArray"]:
        return "SOMASparseNDArray"

    def create(
        self,
        type: pa.DataType,
        shape: Union[NTuple, List[int]],
        platform_config: Optional[PlatformConfig] = None,
    ) -> "SparseNDArray":
        """
        Create a ``SparseNDArray`` named with the URI.

        :param type: an Arrow type defining the type of each element in the array. If the type is unsupported, an error will be raised.

        :param shape: the length of each domain as a list, e.g., [100, 10]. All lengths must be in the positive int64 range.
        """

        # check on shape
        if len(shape) == 0 or any(e <= 0 for e in shape):
            raise ValueError(
                "DenseNDArray shape must be non-zero length tuple of ints > 0"
            )

        if not pa.types.is_primitive(type):
            raise TypeError(
                "Unsupported type - DenseNDArray only supports primtive Arrow types"
            )

        level = self._tiledb_platform_config.string_dim_zstd_level
        create_options = tdbpc.from_param(platform_config).create_options()

        dims = []
        for n, e in enumerate(shape):
            dim_name = f"soma_dim_{n}"
            dim = tiledb.Dim(
                name=dim_name,
                domain=(0, e - 1),
                tile=create_options.dim_tile(dim_name, min(e, 2048)),
                dtype=np.int64,
                filters=create_options.dim_filters(
                    dim_name, [dict(_type="ZstdFilter", level=level)]
                ),
            )
            dims.append(dim)
        dom = tiledb.Domain(dims, ctx=self._ctx)

        attrs = [
            tiledb.Attr(
                name="soma_data",
                dtype=util_arrow.tiledb_type_from_arrow_type(type),
                filters=create_options.attr_filters("soma_data", ["ZstdFilter"]),
                ctx=self._ctx,
            )
        ]

        cell_order, tile_order = create_options.cell_tile_orders()

        # TODO: code-dedupe w/ regard to DenseNDArray. The two creates are
        # almost identical & could share a common parent-class _create() method.
        sch = tiledb.ArraySchema(
            domain=dom,
            attrs=attrs,
            sparse=True,
            allows_duplicates=False,
            offsets_filters=create_options.offsets_filters(),
            capacity=create_options.get("capacity", 100000),
            tile_order=tile_order,
            cell_order=cell_order,
            ctx=self._ctx,
        )

        tiledb.Array.create(self._uri, sch, ctx=self._ctx)

        self._common_create()  # object-type metadata etc

        return self

    @property
    def shape(self) -> NTuple:
        """
        Return length of each dimension, always a list of length ``ndim``
        """
        with self._tiledb_open() as A:
            return cast(NTuple, A.schema.domain.shape)

    def reshape(self, shape: NTuple) -> None:
        raise NotImplementedError("reshape operation not implemented.")

    @property
    def ndim(self) -> int:
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
        return cast(
            int,
            clib.SOMAReader(
                self.uri,
                platform_config={} if self._ctx is None else self._ctx.config().dict(),
            ).nnz(),
        )

    def read_sparse_tensor(
        self,
        coords: SparseNdCoordinates,
        *,
        format: Literal["coo", "csr", "csc"] = "coo",
    ) -> Iterator[Union[pa.SparseCOOTensor, pa.SparseCSCMatrix, pa.SparseCSRMatrix]]:
        """
        Read a use-defined slice of the SparseNDArray and return as an Arrow sparse tensor.

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

        Acceptable ways to index
        ------------------------
        See `read_table`.

        Returns
        -------
        The requested data in the sparse tensor format specified.
        """

        if format != "coo" and self.ndim != 2:
            raise ValueError(f"Format {format} only supported for 2D SparseNDArray")
        if format not in ("coo", "csr", "csc"):
            raise NotImplementedError("format not implemented")

        # The full shape of the data is a template for the shape of the answer.
        # For example if the data is 2D and coords is (None, [5,7,9]), then
        # output shape[0] will be taken from the full shape, and output shape[1] will
        # be 3.
        with self._tiledb_open("r") as A:
            shape = A.shape

        for arrow_tbl in self.read_table(coords):
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

                coo_data = arrow_tbl.column("soma_data").to_numpy()
                coo_coords = np.array(
                    [
                        arrow_tbl.column(f"soma_dim_{n}").to_numpy()
                        for n in range(self.ndim)
                    ]
                ).T
                yield pa.SparseCOOTensor.from_numpy(coo_data, coo_coords, shape=shape)

            elif format in ("csr", "csc"):
                # Temporary: as these must be 2D, convert to scipy COO and use
                # scipy to perform conversions.  C++ reader will be nicer!
                data = arrow_tbl.column("soma_data").to_numpy()
                row = arrow_tbl.column("soma_dim_0").to_numpy()
                col = arrow_tbl.column("soma_dim_1").to_numpy()
                scipy_coo = sp.coo_array((data, (row, col)), shape=shape)
                if format == "csr":
                    yield pa.SparseCSRMatrix.from_scipy(scipy_coo.tocsr())
                if format == "csc":
                    yield pa.SparseCSCMatrix.from_scipy(scipy_coo.tocsc())

    def read_table(self, coords: SparseNdCoordinates) -> Iterator[pa.Table]:
        """
        Read a user-defined slice of the sparse array and return in COO format
        as an Arrow Table

        Acceptable ways to index
        ------------------------

        * None
        * A sequence of coordinates is accepted, one per dimension.
        * Sequence length must be at least one and <= number of dimensions.
        * If the sequence contains missing coordinates (length less than number of dimensions),
          then `slice(None)` -- i.e. no constraint -- is assumed for the missing dimensions.
        * Per-dimension, explicitly specified coordinates can be one of: None, a value, a
          list/ndarray/paarray/etc of values, a slice, etc.
        * Slices are doubly inclusive: slice(2,4) means [2,3,4] not [2,3]. Slice steps can only be +1.
          Slices can be `slice(None)`, meaning select all in that dimension, but may not be half-specified:
          `slice(2,None)` and `slice(None,4)` are both unsupported.
        * Negative indexing is unsupported.
        """

        with self._tiledb_open("r") as A:
            sr = clib.SOMAReader(
                self._uri,
                name=self.__class__.__name__,
                schema=A.schema,
                platform_config={} if self._ctx is None else self._ctx.config().dict(),
            )

            if coords is not None:
                if not isinstance(coords, (list, tuple)):
                    raise TypeError(
                        f"coords type {type(coords)} unsupported; expected list or tuple"
                    )
                if len(coords) < 1 or len(coords) > A.schema.domain.ndim:
                    raise ValueError(
                        f"coords {coords} must have length between 1 and ndim ({A.schema.domain.ndim}); got {len(coords)}"
                    )

                for i, coord in enumerate(coords):
                    #                # Example: coords = [None, 3, slice(4,5)]
                    #                # coord takes on values None, 3, and slice(4,5) in this loop body.
                    dim_name = A.schema.domain.dim(i).name
                    if coord is None:
                        pass  # No constraint; select all in this dimension
                    elif isinstance(coord, int):
                        sr.set_dim_points(dim_name, [coord])
                    elif isinstance(coord, np.ndarray):
                        if coord.ndim != 1:
                            raise ValueError(
                                f"only 1D numpy arrays may be used to index; got {coord.ndim}"
                            )
                        sr.set_dim_points(dim_name, coord)
                    elif isinstance(coord, slice):
                        ned = A.nonempty_domain()  # None iff the array has no data
                        lo_hi = util.slice_to_range(coord, ned[i]) if ned else None
                        if lo_hi is not None:
                            lo, hi = lo_hi
                            if lo < 0 or hi < 0:
                                raise ValueError(
                                    f"slice start and stop may not be negative; got ({lo}, {hi})"
                                )
                            if lo > hi:
                                raise ValueError(
                                    f"slice start must be <= slice stop; got ({lo}, {hi})"
                                )
                            sr.set_dim_ranges(dim_name, [lo_hi])
                        # Else, no constraint in this slot. This is `slice(None)` which is like
                        # Python indexing syntax `[:]`.
                    elif isinstance(
                        coord, (collections.abc.Sequence, pa.Array, pa.ChunkedArray)
                    ):
                        sr.set_dim_points(dim_name, coord)
                    else:
                        raise TypeError(
                            f"coord type {type(coord)} at slot {i} unsupported"
                        )

            sr.submit()

            # This requires careful handling in the no-data case.
            #
            # When there is at least one table which has non-zero length, we could use the following:
            #   while arrow_table := sr.read_next():
            #     yield arrow_table
            # This respects the caller's expectation that there will always be at least one table in
            # the iterator. (For example, pd.concat(self.read_as_pandas(coords)) will raise an
            # exception if the iterator has zero elements.)
            #
            # But when there is only zero-length data available, the above does *not* respect
            # caller expectations: because a zero-length pyarrow.Table is falsy (not truthy)
            # the while-loop becomes zero-pass.
            #
            # A tempting alternative is to instead write:
            #   for arrow_table in sr.read_next():
            #       yield arrow_table
            # This is correctly one-pass. However, it yields an iterator of pyarrow.ChunkedArray,
            # not an iterator of pyarrow.Table.
            #
            # For this reason, we use the following i > 0 check to guarantee a minimum of one
            # pass through the yielded iterator even when the resulting table is zero-length.
            i = 0
            while True:
                arrow_table = sr.read_next()
                if not arrow_table and i > 0:
                    break
                i += 1
                yield arrow_table

    def read_as_pandas(self, coords: SparseNdCoordinates) -> Iterator[pd.DataFrame]:
        """
        Read a user-defined slice of the sparse array and return as a Pandas DataFrame
        containing COO data.
        """
        for arrow_tbl in self.read_table(coords):
            yield arrow_tbl.to_pandas()

    def read_as_pandas_all(
        self, coords: Optional[SparseNdCoordinates] = None
    ) -> pd.DataFrame:
        """
        Return the sparse array as a single Pandas DataFrame containing COO data.
        """
        if coords is None:
            coords = (slice(None),) * self.ndim
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
        Write an Arrow sparse tensor to the SparseNDArray. The coordinates in the Arrow
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
            if self.ndim != 2:
                raise ValueError(
                    f"Unable to write 2D Arrow sparse matrix to {self.ndim}D SparseNDArray"
                )
            # TODO: the `to_scipy` function is not zero copy. Need to explore zero-copy options.
            sp = tensor.to_scipy().tocoo()
            with self._tiledb_open("w") as A:
                A[sp.row, sp.col] = sp.data
            return

        raise TypeError("Unsupported tensor type")

    def write_table(self, arrow_table: pa.Table) -> None:
        """
        Write a COO table, with columns named ``soma_dim_0``, ..., ``soma_dim_N`` and ``soma_data``
        to the dense nD array.
        """
        data = arrow_table.column("soma_data").to_numpy()
        coord_tbl = arrow_table.drop(["soma_data"])
        coords = tuple(
            coord_tbl.column(f"soma_dim_{n}").to_numpy()
            for n in range(coord_tbl.num_columns)
        )
        with self._tiledb_open("w") as A:
            A[coords] = data
