import tiledb
from .soma_options      import SOMAOptions
from .tiledb_group      import TileDBGroup
from .annotation_matrix import AnnotationMatrix
import tiledbsc.util as util

import pandas as pd
import scipy

from typing import Optional, Dict
import os

class AnnotationMatrixGroup(TileDBGroup):
    """
    Nominally for soma obsm and varm.
    """

    dim_name: str


    # ----------------------------------------------------------------
    def __init__(
        self,
        uri: str,
        name: str, # 'obsm' or 'varm'
        parent: Optional[TileDBGroup] = None,
    ):
        """
        See the TileDBObject constructor.
        """
        assert(name in ['obsm', 'varm'])
        super().__init__(uri=uri, name=name, parent=parent)
        self.dim_name = 'obs_id' if name == 'obsm' else 'var_id'

    # ----------------------------------------------------------------
    def from_anndata(self, annotation_matrices, dim_values):
        """
        Populates the obsm/ or varm/ subgroup for a SOMA object, then writes all the components
        arrays under that group.

        :param annotation_matrices: anndata.obsm, anndata.varm, or anndata.raw.varm.
        :param dim_values: anndata.obs_names, anndata.var_names, or anndata.raw.var_names.
        """

        self.open('w')

        for matrix_name in annotation_matrices.keys():
            anndata_matrix = annotation_matrices[matrix_name]
            matrix_uri = os.path.join(self.uri, matrix_name)
            annotation_matrix = AnnotationMatrix(
                uri=matrix_uri,
                name=matrix_name,
                dim_name=self.dim_name,
                parent=self,
            )
            annotation_matrix.from_anndata(anndata_matrix, dim_values)
            self.add(annotation_matrix)
        self.close()

    # ----------------------------------------------------------------
    def to_dict_of_csr(self) -> Dict[str, scipy.sparse.csr_matrix]:
        """
        Reads the obsm/varm group-member arrays into a dict from name to member array.
        Member arrays are returned in sparse CSR format.
        """

        grp = None
        try: # Not all groups have all four of obsm, obsp, varm, and varp.
            grp = tiledb.Group(self.uri, mode='r')
        except:
            pass
        if grp == None:
            if self.verbose:
                print(f"{self.indent}{self.uri} not found")
            return {}

        if self.verbose:
            s = util.get_start_stamp()
            print(f"{self.indent}START  read {self.uri}")

        # TODO: fold this element-enumeration into the TileDB group class.  Maybe on the same PR
        # where we support somagroup['name'] with overloading of the [] operator.
        matrices_in_group = {}
        for element in grp:
            with tiledb.open(element.uri) as A:
                with tiledb.open(element.uri) as A:
                    if self.verbose:
                        s2 = util.get_start_stamp()
                        print(f"{self.indent}START  read {element.uri}")

                    df = pd.DataFrame(A[:])
                    df.set_index(self.dim_name, inplace=True)
                    matrix_name = os.path.basename(element.uri) # e.g. 'X_pca'
                    matrices_in_group[matrix_name] = df.to_numpy()

                    if self.verbose:
                        print(util.format_elapsed(s2, f"{self.indent}FINISH read {element.uri}"))

        grp.close()

        if self.verbose:
            print(util.format_elapsed(s, f"{self.indent}FINISH read {self.uri}"))

        return matrices_in_group
