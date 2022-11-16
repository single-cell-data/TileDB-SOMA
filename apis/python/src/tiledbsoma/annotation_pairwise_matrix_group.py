from typing import Dict, Iterator, Optional, Sequence

import scipy.sparse as sp
import tiledb

import tiledbsoma.util as util

from .annotation_dataframe import AnnotationDataFrame
from .assay_matrix import AssayMatrix
from .logging import log_io
from .tiledb_group import TileDBGroup
from .types import Labels, Matrix


class AnnotationPairwiseMatrixGroup(TileDBGroup):
    """
    Nominally for SOMA ``obsp`` and ``varp``. You can find element names using ``soma.obsp.keys()``;
    you access elements using ``soma.obsp['distances']`` etc., or soma.obsp.distances if you prefer.
    (The latter syntax is possible when the element name doesn't have dashes, dots, etc. in it.)
    """

    def __init__(
        self,
        uri: str,
        name: str,
        row_dataframe: AnnotationDataFrame,  # Nominally a reference to soma.obs
        col_dataframe: AnnotationDataFrame,  # Nominally a reference to soma.var
        *,
        parent: Optional[TileDBGroup] = None,
    ):
        """
        See the ``TileDBObject`` constructor.
        See ``AssayMatrix`` for the rationale behind retaining references to the ``row_dataframe`` and
        ``col_dataframe`` objects.
        """
        assert name in ["obsp", "varp"]
        super().__init__(uri=uri, name=name, parent=parent)
        if name == "obsp":
            self.row_dim_name = "obs_id_i"
            self.col_dim_name = "obs_id_j"
        else:
            self.row_dim_name = "var_id_i"
            self.col_dim_name = "var_id_j"
        self.row_dataframe = row_dataframe
        self.col_dataframe = col_dataframe

    # ----------------------------------------------------------------
    def keys(self) -> Sequence[str]:
        """
        For obsp and varp, ``.keys()`` is a keystroke-saver for the more general group-member
        accessor ``.get_member_names()``.
        """
        return self.get_member_names()

    # ----------------------------------------------------------------
    def __repr__(self) -> str:
        """
        Default display of soma.obsp and soma.varp.
        """
        return ", ".join(f"'{key}'" for key in self.keys())

    # ----------------------------------------------------------------
    def __getattr__(self, name: str) -> Optional[AssayMatrix]:
        """
        This is called on ``soma.obsp.name`` when ``name`` is not already an attribute.
        This way you can do ``soma.obsp.distances`` as an alias for ``soma.obsp['distances']``.
        """
        with self._open() as G:
            if name not in G:
                raise AttributeError(
                    f"'{self.__class__.__name__}' object has no attribute '{name}'"
                )
        return self[name]

    # ----------------------------------------------------------------
    # At the tiledb-py API level, *all* groups are name-indexable.  But here at the tiledbsoma-py
    # level, we implement name-indexing only for some groups:
    #
    # * Most soma member references are done using Python's dot syntax. For example, rather than
    #   soma['X'], we have simply soma.X, and likewise, soma.raw.X.  Likewise soma.obs and soma.var.
    #
    # * Index references are supported for obsm, varm, obsp, varp, and uns. E.g.
    #   soma.obsm['X_pca'] or soma.uns['neighbors']['params']['method']
    #
    # * Overloading the ``[]`` operator at the TileDBGroup level isn't necessary -- e.g. we don't need
    #   soma['X'] when we have soma.X -- but also it causes circular-import issues in Python.
    #
    # * Rather than doing a TileDBIndexableGroup which overloads the ``[]`` operator, we overload
    #   the ``[]`` operator separately in the various classes which need indexing. This is again to
    #   avoid circular-import issues, and means that [] on ``AnnotationMatrixGroup`` will return an
    #   ``AnnotationMatrix, [] on ``UnsGroup`` will return ``UnsArray`` or ``UnsGroup``, etc.
    def __getitem__(self, name: str) -> Optional[AssayMatrix]:
        """
        Returns an ``AssayMatrix`` element at the given name within the group, or ``None`` if no such
        member exists.  Overloads the ``[...]`` operator.
        """

        with self._open("r") as G:
            try:
                obj = G[name]  # This returns a tiledb.object.Object.
            except tiledb.TileDBError:
                return None

            if obj.type == tiledb.tiledb.Group:
                raise Exception(
                    "Internal error: found group element where array element was expected."
                )
            if obj.type != tiledb.libtiledb.Array:
                raise Exception(
                    f"Internal error: found group element neither subgroup nor array: type is {str(obj.type)}"
                )

            return AssayMatrix(
                uri=obj.uri,
                name=name,
                row_dim_name=self.row_dim_name,
                col_dim_name=self.col_dim_name,
                row_dataframe=self.row_dataframe,
                col_dataframe=self.col_dataframe,
                parent=self,
            )

    # ----------------------------------------------------------------
    def __contains__(self, name: str) -> bool:
        """
        Implements ``"namegoeshere" in soma.obsp/soma.varp``.
        """
        with self._open("r") as G:
            return name in G

    # ----------------------------------------------------------------
    def __iter__(self) -> Iterator[AssayMatrix]:
        """
        Implements ``for matrix in soma.obsp: ...`` and ``for matrix in soma.varp: ...``
        """
        for name, uri in self.get_member_names_to_uris().items():
            yield AssayMatrix(
                uri=uri,
                name=name,
                row_dim_name=self.row_dim_name,
                col_dim_name=self.col_dim_name,
                row_dataframe=self.row_dataframe,
                col_dataframe=self.col_dataframe,
                parent=self,
            )

    # ----------------------------------------------------------------
    def remove(self, matrix_name: str) -> None:
        """
        Removes a component of the ``obsp`` or ``varp`` subgroup for a SOMA object.
        Implements ``del soma.obsp['distances']`` etc.
        """
        self._remove_object_by_name(matrix_name)

    def __delattr__(self, matrix_name: str) -> None:
        """
        Removes a component of the ``obsp`` or ``varp`` subgroup for a SOMA object,
        when invoked as ``del soma.obsp.namegoeshere``.
        """
        self.remove(matrix_name)

    def __delitem__(self, matrix_name: str) -> None:
        """
        Removes a component of the ``obsp`` or ``varp`` subgroup for a SOMA object,
        when invoked as ``del soma.obsp["namegoeshere"]``.
        """
        self.remove(matrix_name)

    # ----------------------------------------------------------------
    def add_matrix_from_matrix_and_dim_values(
        self,
        matrix: Matrix,
        dim_values: Labels,
        matrix_name: str,
        *,
        schema_only: bool = False,
    ) -> None:
        """
        Populates a component of the ``obsp`` or ``varp`` subgroup for a SOMA object.

        :param matrix: element of anndata.obsp or anndata.varp.
        :param dim_values: anndata.obs_names or anndata.var_names.
        :param matrix_name_name: name of the matrix, like ``"distances"``.
        """

        # Must be done first, to create the parent directory
        self.create_unless_exists()

        # See comments in that function
        matrix_uri = self._get_child_uri(matrix_name)

        annotation_pairwise_matrix = AssayMatrix(
            uri=matrix_uri,
            name=matrix_name,
            row_dim_name=self.row_dim_name,
            col_dim_name=self.col_dim_name,
            row_dataframe=self.row_dataframe,
            col_dataframe=self.col_dataframe,
            parent=self,
        )
        annotation_pairwise_matrix.from_matrix_and_dim_values(
            matrix,
            dim_values,
            dim_values,
            schema_only=schema_only,
        )
        self._add_object(annotation_pairwise_matrix)

    # ----------------------------------------------------------------
    def to_dict_of_csr(
        self, obs_df_index: Labels, var_df_index: Labels
    ) -> Dict[str, sp.csr_matrix]:
        """
        Reads the ``obsp`` or ``varp`` group-member arrays into a dict from name to member array.
        Member arrays are returned in sparse CSR format.
        """

        grp = None
        try:  # Not all groups have all four of obsm, obsp, varm, and varp.
            grp = tiledb.Group(self.uri, mode="r", ctx=self._ctx)
        except tiledb.TileDBError:
            pass
        if grp is None:
            log_io(None, f"{self._indent}{self.nested_name} not found")
            return {}

        s = util.get_start_stamp()
        log_io(None, f"{self._indent}START  read {self.nested_name}")

        matrices_in_group = {}
        for element in self:
            s2 = util.get_start_stamp()
            log_io(None, f"{self._indent}START  read {element.nested_name}")

            matrix_name = element.name
            matrices_in_group[matrix_name] = element.to_csr_matrix(
                obs_df_index, var_df_index
            )

            log_io(
                None,
                util.format_elapsed(
                    s2, f"{self._indent}FINISH read {element.nested_name}"
                ),
            )

        grp.close()

        log_io(
            f"Read {self.nested_name}",
            util.format_elapsed(s, f"{self._indent}FINISH READING {self.nested_name}"),
        )

        return matrices_in_group
