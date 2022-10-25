from typing import Any, Dict, Literal, Optional, Tuple, Union, cast

import tiledb

from .soma_collection import CollectionBase
from .soma_dataframe import DataFrame
from .soma_dense_nd_array import SOMADenseNdArray
from .soma_indexed_dataframe import IndexedDataFrame
from .soma_sparse_nd_array import SOMASparseNdArray
from .tiledb_object import TileDBObject
from .tiledb_platform_config import TileDBPlatformConfig


class Measurement(CollectionBase[TileDBObject]):
    """
    A ``Measurement`` is a sub-element of a ``Experiment``, and is otherwise a specialized ``Collection`` with pre-defined fields:

    ``var``: ``DataFrame``

    Primary annotations on the variable axis, for variables in this measurement (i.e., annotates columns of ``X``). The contents of the ``soma_rowid`` pseudo-column define the variable index domain, AKA varid. All variables for this measurement must be defined in this dataframe.

    ``X``: ``Collection`` of ``SOMASparseNdArray``

    A collection of sparse matrices, each containing measured feature values. Each matrix is indexed by ``[obsid, varid]``.

    ``obsm``: ``Collection`` of ``SOMADenseNdArray``

    A collection of dense matrices containing annotations of each ``obs`` row. Has the same shape as ``obs``, and is indexed with ``obsid``.

    ``obsp``: ``Collection`` of ``SOMASparseNdArray``

    A collection of sparse matrices containing pairwise annotations of each ``obs`` row. Indexed with ``[obsid_1, obsid_2]``.

    ``varm``: ``Collection`` of ``SOMADenseNdArray``

    A collection of dense matrices containing annotations of each ``var`` row. Has the same shape as ``var``, and is indexed with ``varid``.

    ``varp``: ``Collection`` of ``SOMASparseNdArray``

    A collection of sparse matrices containing pairwise annotations of each ``var`` row. Indexed with ``[varid_1, varid_2]``
    """

    _subclass_constrained_types: Dict[str, Tuple[str, ...]] = {
        "var": ("DataFrame", "IndexedDataFrame"),
        "X": ("Collection",),
        "obsm": ("Collection",),
        "obsp": ("Collection",),
        "varm": ("Collection",),
        "varp": ("Collection",),
    }

    def __init__(
        self,
        uri: str,
        *,
        # Non-top-level objects can have a parent to propagate context, depth, etc.
        parent: Optional[CollectionBase[Any]] = None,
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
    def soma_type(self) -> Literal["Measurement"]:
        return "Measurement"

    def create(self) -> "Measurement":
        """
        Creates the data structure on disk/S3/cloud.
        """
        super().create()
        return self

    @property
    def var(self) -> Union[DataFrame, IndexedDataFrame]:
        return cast(Union[DataFrame, IndexedDataFrame], self["var"])

    @property
    def X(self) -> CollectionBase[Union[SOMADenseNdArray, SOMASparseNdArray]]:
        return cast(
            CollectionBase[Union[SOMADenseNdArray, SOMASparseNdArray]], self["X"]
        )

    @property
    def obsm(self) -> CollectionBase[SOMADenseNdArray]:
        return cast(CollectionBase[SOMADenseNdArray], self["obsm"])

    @property
    def obsp(self) -> CollectionBase[SOMASparseNdArray]:
        return cast(CollectionBase[SOMASparseNdArray], self["obsp"])

    @property
    def varm(self) -> CollectionBase[SOMADenseNdArray]:
        return cast(CollectionBase[SOMADenseNdArray], self["varm"])

    @property
    def varp(self) -> CollectionBase[SOMASparseNdArray]:
        return cast(CollectionBase[SOMASparseNdArray], self["varp"])
