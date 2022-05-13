import tiledb
from .tiledb_array import TileDBArray
from .tiledb_group import TileDBGroup
from .soma_options import SOMAOptions
import tiledbsc.util as util

import pandas

from typing import Optional, Tuple, List

class AnnotationDataFrame(TileDBArray):
    """
    Nominally for obs and var data within a soma. These have one string dimension, and multiple attributes.
    """

    dim_name: str

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
        assert(name in ['obs', 'var'])
        super().__init__(uri=uri, name=name, parent=parent)
        self.dim_name = name + '_id'


    # ----------------------------------------------------------------
    def from_dataframe(self, dataframe: pandas.DataFrame, extent: int) -> None:
        """
        Populates the obs/ or var/ subgroup for a SOMA object.

        :param dataframe: anndata.obs, anndata.var, anndata.raw.var.
        :param extent: TileDB extent parameter for the array schema.
        """

        offsets_filters = tiledb.FilterList(
            [tiledb.PositiveDeltaFilter(), tiledb.ZstdFilter(level=-1)]
        )
        dim_filters = tiledb.FilterList([tiledb.ZstdFilter(level=-1)])
        attr_filters = tiledb.FilterList([tiledb.ZstdFilter(level=-1)])

        if self.verbose:
            s = util.get_start_stamp()
            print(f"{self.indent}START  WRITING {self.uri}")

        # Make the row-names column (barcodes for obs, gene names for var) explicitly named.
        # Otherwise it'll be called '__tiledb_rows'.
        #
        # Before:
        #
        #   >>> anndata.obs
        #                  orig.ident nCount_RNA nFeature_RNA ...
        #   ATGCCAGAACGACT 0          70.0       47           ...
        #   CATGGCCTGTGCAT 0          85.0       52           ...
        #   ...            ...        ...        ...          ...
        #   GGAACACTTCAGAC 0          150.0      30           ...
        #   CTTGATTGATCTTC 0          233.0      76           ...
        #
        # After:
        #
        #   >>> anndata.obs.rename_axis('obs_id')
        #                  orig.ident nCount_RNA nFeature_RNA ...
        #   obs_id
        #   ATGCCAGAACGACT 0          70.0       47           ...
        #   CATGGCCTGTGCAT 0          85.0       52           ...
        #   ...            ...        ...        ...          ...
        #   GGAACACTTCAGAC 0          150.0      30           ...
        #   CTTGATTGATCTTC 0          233.0      76           ...
        dataframe = dataframe.rename_axis(self.dim_name)

        mode = 'ingest'
        if self.exists():
            mode = 'append'
            if self.verbose:
                print(f"{self.indent}Re-using existing array {self.uri}")

        tiledb.from_pandas(
            uri=self.uri,
            dataframe=dataframe,
            name=self.name,
            sparse=True,
            allows_duplicates=False,
            offsets_filters=offsets_filters,
            attr_filters=attr_filters,
            dim_filters=dim_filters,
            capacity=100000,
            tile=extent,
            ctx=self.ctx,
            mode=mode,
        )

        if self.verbose:
            print(util.format_elapsed(s, f"{self.indent}FINISH WRITING {self.uri}"))


    # ----------------------------------------------------------------
    def to_dataframe(self) -> pandas.DataFrame:
        """
        Reads the TileDB obs or var array and returns a type of pandas dataframe
        and dimension values.
        """

        if self.verbose:
            s = util.get_start_stamp()
            print(f"{self.indent}START  read {self.uri}")

        with tiledb.open(self.uri) as A:
            # We could use A.df[:] to set the index_name to 'obs_id' or 'var_id'.
            # However, the resulting dataframe has obs_id/var_id as strings, not
            # bytes, resulting in `KeyError` elsewhere in the code.
            df = pandas.DataFrame(A[:])
            df = df.set_index(self.dim_name)

        if self.verbose:
            print(util.format_elapsed(s, f"{self.indent}FINISH read {self.uri}"))

        return df

    # ----------------------------------------------------------------
    def get_dim_values(self) -> List[str]:
        """
        TODO
        """
        with tiledb.open(self.uri) as A:
            df = A[:]
            return df[self.dim_name].tolist()
