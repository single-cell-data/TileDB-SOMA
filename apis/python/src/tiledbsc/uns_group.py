import tiledb

from .soma_options import SOMAOptions
from .tiledb_group import TileDBGroup
from .uns_array    import UnsArray
import tiledbsc.util as util

import anndata as ad
import pandas  as pd
import scipy
import numpy

from typing import Optional, Dict

import os

class UnsGroup(TileDBGroup):
    """
    Nominally for soma uns.
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
    def from_anndata_uns(self, uns: ad.compat.OverloadedDict):
        """
        Populates the uns group for the soma object.

        :param uns: anndata.uns.
        """

        if self.verbose:
            s = util.get_start_stamp()
            print(f"{self.indent}START  WRITING {self.uri}")

        self.open("w")

        for key in uns.keys():
            component_uri = os.path.join(self.uri, key)
            value = uns[key]

            if key == 'rank_genes_groups':
                # TODO:
                # This is of type 'structured array':
                # https://numpy.org/doc/stable/user/basics.rec.html
                #
                # >>> a.uns['rank_genes_groups']['names'].dtype
                # dtype([('0', 'O'), ('1', 'O'), ('2', 'O'), ('3', 'O'), ('4', 'O'), ('5', 'O'), ('6', 'O'), ('7', 'O')])
                # >>> type(a.uns['rank_genes_groups']['names'])
                # <class 'numpy.ndarray'>
                #
                # We don’t have a way to model this directly in TileDB schema right now. We support
                # multiplicities of a single scalar type, e.g. a record array with cell_val_num==3
                # and float32 slots (which would correspond to numpy record array
                # np.dtype([("field1", "f4"), ("field2", "f4"), ("field3", "f4",)])). We don’t
                # support nested cells, AKA "list" type.
                #
                # This could, however, be converted to a dataframe and ingested that way.
                print("      Skipping structured array:", component_uri)
                continue

            if isinstance(value, dict) or isinstance(value, ad.compat.OverloadedDict):
                # Nested data, e.g. a.uns['draw-graph']['params']['layout']
                subgroup = UnsGroup(uri=component_uri, name=key, parent=self)
                subgroup.from_anndata_uns(value)
                self.add(subgroup)
                continue

            array = UnsArray(uri=component_uri, name=key, parent=self)

            if isinstance(value, pd.DataFrame):
                array.from_pandas_dataframe(value)
                self.add(array)

            elif isinstance(value, scipy.sparse.csr_matrix):
                array.from_scipy_csr(value)
                self.add(array)

            elif array.maybe_from_numpyable_object(value):
                self.add(array)

            else:
                print("      Skipping unrecognized type:", component_uri, type(value))

        self.close()

        if self.verbose:
            print(util.format_elapsed(s, f"{self.indent}FINISH WRITING {self.uri}"))

    # ----------------------------------------------------------------
    def to_dict_of_matrices(self) -> Dict:
        """
        Reads the recursive group/array uns data from TileDB storage and returns them as a recursive dict of matrices.
        """
        grp = None
        try: # Not all groups have uns
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
        retval = {}
        for element in grp:
            name = os.path.basename(element.uri) # TODO: update for tiledb cloud

            if element.type == tiledb.tiledb.Group:
                child_group = UnsGroup(uri=element.uri, name=name, parent=self)
                retval[name] = child_group.to_dict_of_matrices()

            elif element.type == tiledb.libtiledb.Array:
                child_array = UnsArray(uri=element.uri, name=name, parent=self)
                retval[name] = child_array.to_matrix()

            else:
                raise Exception("Internal error: found uns group element neither group nor array: type is", str(element.type))

        grp.close()

        if self.verbose:
            print(util.format_elapsed(s, f"{self.indent}FINISH read {self.uri}"))

        return retval
