import os
from typing import Dict, List, Optional

import pandas as pd
import scipy.sparse
import tiledb

import tiledbsc.util as util

from .annotation_matrix import AnnotationMatrix
from .logging import logger
from .tiledb_group import TileDBGroup


class AnnotationMatrixGroup(TileDBGroup):
    """
    Nominally for soma obsm and varm. You can find element names using soma.obsm.keys(); you access
    elements using soma.obsm['X_pca'] etc., or soma.obsm.X_pca if you prefer.  (The latter syntax is
    possible when the element name doesn't have dashes, dots, etc. in it.)
    """

    dim_name: str

    # ----------------------------------------------------------------
    def __init__(
        self,
        uri: str,
        name: str,  # 'obsm' or 'varm'
        parent: Optional[TileDBGroup] = None,
    ):
        """
        See the TileDBObject constructor.
        """
        assert name in ["obsm", "varm"]
        super().__init__(uri=uri, name=name, parent=parent)
        self.dim_name = "obs_id" if name == "obsm" else "var_id"

    # ----------------------------------------------------------------
    def keys(self) -> List[str]:
        """
        For `obsm` and `varm`, `.keys()` is a keystroke-saver for the more general group-member
        accessor `._get_member_names()`.
        """
        return self._get_member_names()

    # ----------------------------------------------------------------
    def __repr__(self) -> str:
        """
        Default display of soma.obsm and soma.varm.
        """
        return ", ".join(f"'{key}'" for key in self.keys())

    # ----------------------------------------------------------------
    def __iter__(self) -> List[AnnotationMatrix]:
        """
        Implements `for matrix in soma.obsm: ...` and `for matrix in soma.varm: ...`
        """
        retval = []
        for name, uri in self._get_member_names_to_uris().items():
            matrix = AnnotationMatrix(
                uri=uri, name=name, dim_name=self.dim_name, parent=self
            )
            retval.append(matrix)
        return iter(retval)

    # ----------------------------------------------------------------
    def __getattr__(self, name) -> AnnotationMatrix:
        """
        This is called on `soma.obsm.name` when `name` is not already an attribute.
        This way you can do `soma.obsm.X_tsne` as an alias for `soma.obsm['X_tsne']`.
        """
        with self._open() as G:
            if name not in G:
                raise AttributeError(
                    f"'{self.__class__.__name__}' object has no attribute '{name}'"
                )
        return self[name]

    # ----------------------------------------------------------------
    def add_matrix_from_matrix_and_dim_values(
        self,
        matrix,
        dim_values,
        matrix_name: str,
    ) -> None:
        """
        Populates a component of the `obsm` or `varm` subgroup for a SOMA object.

        :param matrix: element of anndata.obsm, anndata.varm, or anndata.raw.varm.
        :param dim_values: anndata.obs_names, anndata.var_names, or anndata.raw.var_names.
        :param matrix_name_name: name of the matrix, like `"X_tsne"` or `"PCs"`.
        """

        # Must be done first, to create the parent directory
        self.create_unless_exists()

        # See comments in that function
        matrix_uri = self._get_child_uri(matrix_name)

        annotation_matrix = AnnotationMatrix(
            uri=matrix_uri,
            name=matrix_name,
            dim_name=self.dim_name,
            parent=self,
        )
        annotation_matrix.from_matrix_and_dim_values(matrix, dim_values)
        self._add_object(annotation_matrix)

    # ----------------------------------------------------------------
    def to_dict_of_csr(self) -> Dict[str, scipy.sparse.csr_matrix]:
        """
        Reads the obsm/varm group-member arrays into a dict from name to member array.
        Member arrays are returned in sparse CSR format.
        """

        if (
            not self.exists()
        ):  # Not all groups have all four of obsm, obsp, varm, and varp.
            logger.info(f"{self._indent}{self.uri} not found")
            return {}

        s = util.get_start_stamp()
        logger.info(f"{self._indent}START  read {self.uri}")

        with self._open() as G:
            matrices_in_group = {}
            for element in G:
                s2 = util.get_start_stamp()
                logger.info(f"{self._indent}START  read {element.uri}")

                with tiledb.open(element.uri, ctx=self._ctx) as A:
                    df = pd.DataFrame(A[:])
                    df.set_index(self.dim_name, inplace=True)
                    matrix_name = os.path.basename(element.uri)  # e.g. 'X_pca'
                    matrices_in_group[matrix_name] = df.to_numpy()

                logger.info(
                    util.format_elapsed(s2, f"{self._indent}FINISH read {element.uri}")
                )

        logger.info(util.format_elapsed(s, f"{self._indent}FINISH read {self.uri}"))

        return matrices_in_group

    # At the tiledb-py API level, *all* groups are name-indexable.  But here at the tiledbsc-py
    # level, we implement name-indexing only for some groups:
    #
    # * Most soma member references are done using Python's dot syntax. For example, rather than
    #   soma['X'], we have simply soma.X, and likewise, soma.raw.X.  Likewise soma.obs and soma.var.
    #
    # * Index references are supported for obsm, varm, obsp, varp, and uns. E.g.
    #   soma.obsm['X_pca'] or soma.uns['neighbors']['params']['method']
    #
    # * Overloading the `[]` operator at the TileDBGroup level isn't necessary -- e.g. we don't need
    #   soma['X'] when we have soma.X -- but also it causes circular-import issues in Python.
    #
    # * Rather than doing a TileDBIndexableGroup which overloads the `[]` operator, we overload
    #   the `[]` operator separately in the various classes which need indexing. This is again to
    #   avoid circular-import issues, and means that [] on `AnnotationMatrixGroup` will return an
    #   `AnnotationMatrix, [] on `UnsGroup` will return `UnsArray` or `UnsGroup`, etc.
    def __getitem__(self, name) -> AnnotationMatrix:
        """
        Returns an `AnnotationMatrix` element at the given name within the group, or None if no such
        member exists.  Overloads the `[...]` operator.
        """

        with self._open("r") as G:
            if name not in G:
                return None
            obj = G[name]  # This returns a tiledb.object.Object.
            if obj.type == tiledb.tiledb.Group:
                raise Exception(
                    "Internal error: found group element where array element was expected."
                )
            if obj.type != tiledb.libtiledb.Array:
                raise Exception(
                    f"Internal error: found group element neither subgroup nor array: type is {str(obj.type)}"
                )
            return AnnotationMatrix(
                uri=obj.uri, name=name, dim_name=self.dim_name, parent=self
            )

    def __contains__(self, name) -> bool:
        """
        Implements the `in` operator, e.g. `"namegoeshere" in soma.obsm/soma.varm`.
        """
        with self._open("r") as G:
            return name in G
