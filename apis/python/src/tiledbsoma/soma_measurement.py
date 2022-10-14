from typing import Any, Dict, Literal, Optional, Tuple, Union, cast

import tiledb

from .soma_collection import SOMACollectionBase
from .soma_dataframe import SOMADataFrame
from .soma_dense_nd_array import SOMADenseNdArray
from .soma_indexed_dataframe import SOMAIndexedDataFrame
from .soma_sparse_nd_array import SOMASparseNdArray
from .tiledb_object import TileDBObject
from .tiledb_platform_config import TileDBPlatformConfig


class SOMAMeasurement(SOMACollectionBase[TileDBObject]):
    """
    A ``SOMAMeasurement`` is a sub-element of a ``SOMAExperiment``, and is otherwise a specialized ``SOMACollection`` with pre-defined fields:

    ``var``: ``SOMADataFrame``

    Primary annotations on the variable axis, for variables in this measurement (i.e., annotates columns of ``X``). The contents of the ``soma_rowid`` pseudo-column define the variable index domain, AKA varid. All variables for this measurement must be defined in this dataframe.

    ``X``: ``SOMACollection`` of ``SOMASparseNdArray``

    A collection of sparse matrices, each containing measured feature values. Each matrix is indexed by ``[obsid, varid]``.

    ``obsm``: ``SOMACollection`` of ``SOMADenseNdArray``

    A collection of dense matrices containing annotations of each ``obs`` row. Has the same shape as ``obs``, and is indexed with ``obsid``.

    ``obsp``: ``SOMACollection`` of ``SOMASparseNdArray``

    A collection of sparse matrices containing pairwise annotations of each ``obs`` row. Indexed with ``[obsid_1, obsid_2]``.

    ``varm``: ``SOMACollection`` of ``SOMADenseNdArray``

    A collection of dense matrices containing annotations of each ``var`` row. Has the same shape as ``var``, and is indexed with ``varid``.

    ``varp``: ``SOMACollection`` of ``SOMASparseNdArray``

    A collection of sparse matrices containing pairwise annotations of each ``var`` row. Indexed with ``[varid_1, varid_2]``
    """

    _subclass_constrained_types: Dict[str, Tuple[str, ...]] = {
        "var": ("SOMADataFrame", "SOMAIndexedDataFrame"),
        "X": ("SOMACollection",),
        "obsm": ("SOMACollection",),
        "obsp": ("SOMACollection",),
        "varm": ("SOMACollection",),
        "varp": ("SOMACollection",),
    }

    def __init__(
        self,
        uri: str,
        *,
        # Non-top-level objects can have a parent to propagate context, depth, etc.
        parent: Optional[SOMACollectionBase[Any]] = None,
        # Top-level objects should specify these:
        tiledb_platform_config: Optional[TileDBPlatformConfig] = None,
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
    def soma_type(self) -> Literal["SOMAMeasurement"]:
        return "SOMAMeasurement"

    def create(self) -> "SOMAMeasurement":
        """
        Creates the data structure on disk/S3/cloud.
        """
        super().create()
        return self

    @property
    def var(self) -> Union[SOMADataFrame, SOMAIndexedDataFrame]:
        return cast(Union[SOMADataFrame, SOMAIndexedDataFrame], self["var"])

    @property
    def X(self) -> SOMACollectionBase[Union[SOMADenseNdArray, SOMASparseNdArray]]:
        return cast(
            SOMACollectionBase[Union[SOMADenseNdArray, SOMASparseNdArray]], self["X"]
        )

    @property
    def obsm(self) -> SOMACollectionBase[SOMADenseNdArray]:
        return cast(SOMACollectionBase[SOMADenseNdArray], self["obsm"])

    @property
    def obsp(self) -> SOMACollectionBase[SOMASparseNdArray]:
        return cast(SOMACollectionBase[SOMASparseNdArray], self["obsp"])

    @property
    def varm(self) -> SOMACollectionBase[SOMADenseNdArray]:
        return cast(SOMACollectionBase[SOMADenseNdArray], self["varm"])

    @property
    def varp(self) -> SOMACollectionBase[SOMASparseNdArray]:
        return cast(SOMACollectionBase[SOMASparseNdArray], self["varp"])
