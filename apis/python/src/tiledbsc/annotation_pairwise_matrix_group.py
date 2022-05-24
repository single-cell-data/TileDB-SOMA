import tiledb
from .soma_options import SOMAOptions
from .tiledb_group import TileDBGroup
from .annotation_pairwise_matrix import AnnotationPairwiseMatrix
from .assay_matrix import AssayMatrix
import tiledbsc.util as util

import pandas as pd
import scipy

from typing import Optional, Dict, List
import os


class AnnotationPairwiseMatrixGroup(TileDBGroup):
    """
    Nominally for soma obsp and varp.
    """

    row_dim_name: str
    col_dim_name: str

    def __init__(
        self,
        uri: str,
        name: str,
        parent: Optional[TileDBGroup] = None,
    ):
        """
        See the TileDBObject constructor.
        """
        assert name in ["obsp", "varp"]
        super().__init__(uri=uri, name=name, parent=parent)
        if name == "obsp":
            self.row_dim_name = "obs_id_i"
            self.col_dim_name = "obs_id_j"
        else:
            self.row_dim_name = "var_id_i"
            self.col_dim_name = "var_id_j"

    # ----------------------------------------------------------------
    def keys(self):
        """
        For obsp and varp, `.keys()` is a keystroke-saver for the more general group-member
        accessor `._get_member_names()`.
        """
        return self._get_member_names()

    # ----------------------------------------------------------------
    def __iter__(self) -> List[AnnotationPairwiseMatrix]:
        """
        Implements 'for matrix in soma.obsp: ...' and 'for matrix in soma.varp: ...'
        """
        retval = []
        for name, uri in self._get_member_names_to_uris().items():
            matrix = AnnotationPairwiseMatrix(
                uri=uri,
                name=name,
                row_dim_name=self.row_dim_name,
                col_dim_name=self.col_dim_name,
                parent=self,
            )
            retval.append(matrix)
        return iter(retval)

    # ----------------------------------------------------------------
    def from_anndata(self, annotation_pairwise_matrices, dim_values):
        """
        Populates the obsp/ or varp/ subgroup for a SOMA object, then writes all the components
        arrays under that group.

        :param annotation_pairwise_matrices: anndata.obsp, anndata.varp, or anndata.raw.varp.
        :param dim_values: anndata.obs_names, anndata.var_names, or anndata.raw.var_names.
        """

        # Must be done first, to create the parent directory
        self._create()
        for matrix_name in annotation_pairwise_matrices.keys():
            anndata_matrix = annotation_pairwise_matrices[matrix_name]
            matrix_uri = os.path.join(self.uri, matrix_name)
            annotation_pairwise_matrix = AssayMatrix(
                uri=matrix_uri,
                name=matrix_name,
                row_dim_name=self.row_dim_name,
                col_dim_name=self.col_dim_name,
                parent=self,
            )
            annotation_pairwise_matrix.from_matrix(
                anndata_matrix, dim_values, dim_values
            )
            self._add_object(annotation_pairwise_matrix)

    # ----------------------------------------------------------------
    def to_dict_of_csr(self) -> Dict[str, scipy.sparse.csr_matrix]:
        """
        Reads the obsm/varm group-member arrays into a dict from name to member array.
        Member arrays are returned in sparse CSR format.
        """

        grp = None
        try:  # Not all groups have all four of obsm, obsp, varm, and varp.
            grp = tiledb.Group(self.uri, mode="r")
        except:
            pass
        if grp == None:
            if self._verbose:
                print(f"{self._indent}{self.uri} not found")
            return {}

        if self._verbose:
            s = util.get_start_stamp()
            print(f"{self._indent}START  read {self.uri}")

        # TODO: fold this element-enumeration into the TileDB group class.  Maybe on the same PR
        # where we support somagroup['name'] with overloading of the [] operator.
        matrices_in_group = {}
        for element in grp:
            with tiledb.open(element.uri) as A:
                with tiledb.open(element.uri) as A:
                    if self._verbose:
                        s2 = util.get_start_stamp()
                        print(f"{self._indent}START  read {element.uri}")

                    df = pd.DataFrame(A[:])
                    matrix_name = os.path.basename(
                        element.uri
                    )  # TODO: fix for tiledb cloud
                    matrices_in_group[matrix_name] = scipy.sparse.coo_matrix(df).tocsr()
                    # TODO: not working yet:
                    # TypeError: no supported conversion for types: (dtype('O'),)

                    if self._verbose:
                        print(
                            util.format_elapsed(
                                s2, f"{self._indent}FINISH read {element.uri}"
                            )
                        )

        grp.close()

        if self._verbose:
            print(util.format_elapsed(s, f"{self._indent}FINISH read {self.uri}"))

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
    def __getitem__(self, name):
        """
        Returns an `AnnotationPairwiseMatrix` element at the given name within the group, or None if no such
        member exists.  Overloads the [...] operator.
        """

        with self._open("r") as G:
            try:
                obj = G[name]  # This returns a tiledb.object.Object.
            except:
                return None

            if obj.type == tiledb.tiledb.Group:
                raise Exception(
                    "Internal error: found group element where array element was expected."
                )
            if obj.type != tiledb.libtiledb.Array:
                raise Exception(
                    f"Internal error: found group element neither subgroup nor array: type is {str(obj.type)}"
                )

            return AnnotationPairwiseMatrix(
                uri=obj.uri,
                name=name,
                row_dim_name=self.row_dim_name,
                col_dim_name=self.col_dim_name,
                parent=self,
            )

    def __contains__(self, name):
        """
        Implements '"namegoeshere" in soma.obsp/soma.varp'.
        """
        # TODO: this will get easier once TileDB.group.Group supports `name` in `__contains__`.
        # See SC-18057 and https://github.com/single-cell-data/TileDB-SingleCell/issues/113.
        with self._open("r") as G:
            answer = False
            try:
                # This returns a tiledb.object.Object.
                G[name]
                return True
            except:
                return False
