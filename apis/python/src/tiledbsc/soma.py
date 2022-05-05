import os
from typing import Optional, Union

import anndata as ad
import numpy   as np
import pandas  as pd
import pyarrow as pa
import scanpy
import scipy
import tiledb
import tiledbsc.util     as util

class SOMA():
    """ Single-cell group
    A simple class for ingestion of anndata to a TileDB group.

    * The write path (anndata/10X -> TileDB) is exercised here
    * A read path (TileDB -> anndata/10X) is a WIP.
    * The class structure in this source file is likely to change but the
      TileDB Groups this code creates should remain relatively stable, modulo
      updates to the Matrix API itself and/or any errors encountered during testing.
    * See also desc-ann.py in this directory for helpful information to
      reveal the diversity/variety of HDF5 files we process.
    """

    uri: str
    verbose: bool
    config: tiledb.Config
    ctx: tiledb.Ctx

    write_X_chunked_if_csr: bool
    goal_chunk_nnz: int

    # ----------------------------------------------------------------
    def __init__(self, uri: str, verbose: bool = True, config: Optional[tiledb.Config] = None, ctx: Optional[tiledb.Ctx] = None):
        """
        @description Create a new SOMA object. The existing array group is
          opened at the specified array `uri` if one is present, otherwise a new
          array group is created.
        @param uri URI of the TileDB group
        @param verbose Print status messages
        """

        self.uri = uri
        self.verbose = verbose
        self.config = config
        self.ctx = ctx

        if self.ctx is None and self.config is not None:
            self.ctx = tiledb.Ctx(self.config)

        # TODO: make a user-accessible default/override setup
        # See also https://github.com/single-cell-data/TileDB-SingleCell/issues/27
        self.write_X_chunked_if_csr = True
        self.goal_chunk_nnz = 10000000

        # If URI is "/something/test1" then:
        # * obs_uri  is "/something/test1/obs"
        # * var_uri  is "/something/test1/var"
        # * data_uri is "/something/test1/X"

        # If URI is "tiledb://namespace/s3://bucketname/something/test1" then:
        # * obs_uri  is "tiledb://namespace/s3://bucketname/something/test1/obs"
        # * var_uri  is "tiledb://namespace/s3://bucketname/something/test1/var"
        # * data_uri is "tiledb://namespace/s3://bucketname/something/test1/X"

    def set_write_X_chunked_if_csr(self, val: bool):
        """
        Allows the user to disable the default setting which is that X matrices in CSR format
        are written chunked/fragmented as a memory-reduction strategy. If this is set to False,
        X matrices will be converted to COO and ingested all at one go.
        """
        self.write_X_chunked_if_csr = val
        return self

    # ----------------------------------------------------------------
    def from_anndata(self, anndata: ad.AnnData):
        """
        Factory function to instantiate a SOMA object from an input anndata.AnnData object.
        """
        if self.verbose:
            s = util.get_start_stamp()
            print(f"START  SOMA.from_ann")

        anndata = self.decategoricalize(anndata)
        self.write_tiledb_group(anndata)

        if self.verbose:
            print(util.format_elapsed(s,"FINISH  SOMA.from_ann"))

    # ----------------------------------------------------------------
    def from_h5ad(self, input_path: str):
        """
        Factory function to instantiate a SOMA object from an input .h5ad file.
        """
        if self.verbose:
            s = util.get_start_stamp()
            print(f"START  SOMA.from_h5ad {input_path} -> {self.uri}")

        anndata = self.read_h5ad(input_path)
        anndata = self.decategoricalize(anndata)
        self.write_tiledb_group(anndata)

        if self.verbose:
            print(util.format_elapsed(s, f"FINISH SOMA.from_h5ad {input_path} -> {self.uri}"))

    # ----------------------------------------------------------------
    def from_10x(self, input_path: str):
        """
        Factory function to instantiate a SOMA object from an input 10X file.
        """
        if self.verbose:
            s = util.get_start_stamp()
            print(f"START  SOMA.from_10x {input_path} -> {self.uri}")

        anndata = self.read_10x(input_path)

        anndata = self.decategoricalize(anndata)
        self.write_tiledb_group(anndata)

        if self.verbose:
            print(util.format_elapsed(s, f"FINISH SOMA.from_10x {input_path} -> {self.uri}"))

    # ----------------------------------------------------------------
    def read_h5ad(self, input_path: str):
        """
        File-ingestor for .h5ad files
        """
        if self.verbose:
            s = util.get_start_stamp()
            print(f"  START  READING {input_path}")
        anndata = ad.read_h5ad(input_path)
        anndata.obs_names_make_unique()
        anndata.var_names_make_unique()
        if self.verbose:
            print(util.format_elapsed(s, f"  FINISH READING {input_path}"))
        return anndata

    # ----------------------------------------------------------------
    def read_10x(self, input_path: str):
        """
        File-ingestor for 10X files
        """
        if self.verbose:
            s = util.get_start_stamp()
            print(f"  START  READING {input_path}")
        anndata = scanpy.read_10x_h5(input_path)
        anndata.obs_names_make_unique()
        anndata.var_names_make_unique()
        if self.verbose:
            print(util.format_elapsed(s, f"  FINISH READING {input_path}"))
        return anndata

    # ----------------------------------------------------------------
    def decategoricalize(self, anndata: ad.AnnData):
        """
        Performs an in-place typecast for the categorical datatype that pandas can handle.
        Categorical strings -> string; bool -> uint8.
        The uns dataset is deferred.
        """

        if self.verbose:
            s = util.get_start_stamp()
            print(f"  START  DECATEGORICALIZING")

        new_obs = pd.DataFrame.from_dict({k: util.decategoricalize_array(v) for k, v in anndata.obs.items()})
        new_var = pd.DataFrame.from_dict({k: util.decategoricalize_array(v) for k, v in anndata.var.items()})
        for key in anndata.obsm.keys():
            anndata.obsm[key] = util.decategoricalize_array(anndata.obsm[key])
        for key in anndata.varm.keys():
            anndata.varm[key] = util.decategoricalize_array(anndata.varm[key])
        for key in anndata.obsp.keys():
            anndata.obsp[key] = util.decategoricalize_array(anndata.obsp[key])
        for key in anndata.varp.keys():
            anndata.varp[key] = util.decategoricalize_array(anndata.varp[key])

        if anndata.raw == None: # Some datasets have no raw.
            newraw = None
        else:
            # Note there is some code-duplication here between cooked & raw.  However anndata.raw
            # has var not directly assignable ('AttributeError: can't set attribute'), and
            # anndata.AnnData and anndata.Raw have different constructor syntaxes, and raw doesn't
            # have obs or obsm or obsp -- so, it turns out to be simpler to just repeat ourselves a
            # little.

            newvar = pd.DataFrame.from_dict({k: util.decategoricalize_array(v) for k, v in anndata.raw.var.items()})
            for key in anndata.raw.varm.keys():
                anndata.raw.varm[key] = util.decategoricalize_array(anndata.raw.varm[key])
            newraw = ad.Raw(
                anndata,
                X=anndata.raw.X,
                var=newvar,
                varm=anndata.raw.varm,
            )

        anndata = ad.AnnData(
            X=anndata.X,
            dtype=anndata.X.dtype,
            obs=new_obs,
            var=new_var,
            obsm=anndata.obsm,
            obsp=anndata.obsp,
            varm=anndata.varm,
            varp=anndata.varp,
            raw=newraw,
        )

        if self.verbose:
            print(util.format_elapsed(s, f"  FINISH DECATEGORICALIZING"))

        return anndata

    # ----------------------------------------------------------------
    # Intended structure:
    #
    # soma: group
    # |
    # +-- X: group
    # |   +-- data: array
    # |
    # +-- obs: array
    # |
    # +-- var: array
    # |
    # +-- obsm: group
    # |   +-- omfoo: array
    # |   +-- ombar: array
    # |
    # +-- varm: group
    # |   +-- vmfoo: array
    # |   +-- vmbar: array
    # |
    # +-- obsp: group
    # |   +-- opfoo: array
    # |   +-- opbar: array
    # |
    # +-- varp: group
    # |   +-- vpfoo: array
    # |   +-- vpbar: array
    # |
    # +-- raw: group
    #     |
    #     +-- X: group
    #     |   +-- data: array
    #     |
    #     +-- var: array
    #     |
    #     +-- varm: group
    #     |   +-- vmfoo: array
    #     |   +-- vmbar: array
    #     |
    #     +-- varp: group
    #         +-- vpfoo: array
    #         +-- vpbar: array

    def write_tiledb_group(self, anndata: ad.AnnData):
        """
        Top-level writer method for creating a TileDB group for a SOMA object.
        """
        if self.verbose:
            s = util.get_start_stamp()
            print(f"  START  WRITING {self.uri}")

        # ----------------------------------------------------------------
        # Must be done first, to create the parent directory
        tiledb.group_create(self.uri, ctx=self.ctx)
        base_group = tiledb.Group(self.uri, mode="w", ctx=self.ctx)

        # ----------------------------------------------------------------
        X_uri = self.write_X(anndata)
        base_group.add(uri=X_uri, relative=False, name="X")

        # ----------------------------------------------------------------
        obs_uri = self.write_obs_or_var(anndata.obs, "obs", 256)
        base_group.add(uri=obs_uri, relative=False, name="obs")

        var_uri = self.write_obs_or_var(anndata.var, "var", 2048)
        base_group.add(uri=var_uri, relative=False, name="var")

        # ----------------------------------------------------------------
        if len(anndata.obsm.keys()) > 0:
            obsm_uri = self.write_annotation_matrices(anndata.obsm, "obsm", "obs_id", anndata.obs_names)
            base_group.add(uri=obsm_uri, relative=False, name="obsm")

        if len(anndata.varm.keys()) > 0:
            varm_uri = self.write_annotation_matrices(anndata.varm, "varm", "var_id", anndata.var_names)
            base_group.add(uri=varm_uri, relative=False, name="varm")

        if len(anndata.obsp.keys()) > 0:
            obsp_uri = self.write_annotation_pairwise_matrices(anndata.obsp, "obsp", "obs_id", anndata.obs_names)
            base_group.add(uri=obsp_uri, relative=False, name="obsp")

        if len(anndata.varp.keys()) > 0:
            varp_uri = self.write_annotation_pairwise_matrices(anndata.varp, "varp", "var_id", anndata.var_names)
            base_group.add(uri=varp_uri, relative=False, name="varp")

        # ----------------------------------------------------------------
        if self.verbose:
            print(util.format_elapsed(s, f"  FINISH WRITING {self.uri}"))

        base_group.close()

    # ----------------------------------------------------------------
    def write_X(self, anndata: ad.AnnData):
        """
        Populates the X/ subgroup for a SOMA object.
        """
        X_uri = os.path.join(self.uri, "X")
        tiledb.group_create(X_uri, ctx=self.ctx)
        X_group = tiledb.Group(X_uri, mode="w", ctx=self.ctx)

        X_data_uri = self.write_X_array(anndata.X, X_uri, "data", anndata.obs.index, anndata.var.index)
        X_group.add(uri=X_data_uri, relative=False, name="data")

        has_raw = False
        try:
            anndata.raw.X.shape
            has_raw = True
        except:
            pass
        if has_raw:
            X_raw_uri = self.write_X_array(anndata.raw.X, X_uri, "raw", anndata.raw.obs_names, anndata.raw.var_names)
            X_group.add(uri=X_raw_uri, relative=False, name="raw")

        X_group.close()

        return X_uri

    # ----------------------------------------------------------------
    def write_X_array(self, x, group_uri:str, arrayname: str, obs_names, var_names):
        """
        Populates the X/data or X/raw array.
        :param x: is anndata.X or raw
        :param group_uri: is the URI of the parent group, e.g. 'foo/X'
        :param arrayname: is the name of the array within the parent group, e.g. 'data' for 'foo/X/data'
        :param obs_names: and var_names are the names for the axes
        """
        X_array_uri = os.path.join(group_uri, arrayname)
        if self.verbose:
            s = util.get_start_stamp()
            print(f"    START  WRITING {X_array_uri} from {type(x)}")

        self.__create_coo_array(uri=X_array_uri, dim_labels=["obs_id", "var_id"], attr_name="value", mat_dtype=x.dtype)

        # TODO: add chunked support for CSC
        if isinstance(x, scipy.sparse._csr.csr_matrix) and self.write_X_chunked_if_csr:
            self.__ingest_coo_data_rows_chunked(X_array_uri, x, obs_names, var_names)
        else:
            self.__ingest_coo_data_whole(X_array_uri, x, obs_names, var_names)

        if self.verbose:
            print(util.format_elapsed(s, f"    FINISH WRITING {X_array_uri}"))
        return X_array_uri

    # ----------------------------------------------------------------
    def write_obs_or_var(self, obs_or_var_data, obs_or_var_name: str, extent: int):
        """
        Populates the obs/ or var/ subgroup for a SOMA object.
        First argument is anndata.obs or anndata.var; second is "obs" or "var".  In the reference
        pbmc3k_processed dataset, these are of type pandas.core.frame.DataFrame. In further
        testing we may need to switch on the datatype.
        """

        offsets_filters = tiledb.FilterList(
            [tiledb.PositiveDeltaFilter(), tiledb.ZstdFilter(level=-1)]
        )
        dim_filters = tiledb.FilterList([tiledb.ZstdFilter(level=-1)])
        attr_filters = tiledb.FilterList([tiledb.ZstdFilter(level=-1)])

        obs_or_var_uri = os.path.join(self.uri, obs_or_var_name)
        if self.verbose:
            s = util.get_start_stamp()
            print(f"    START  WRITING {obs_or_var_uri}")

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
        obs_or_var_data = obs_or_var_data.rename_axis(obs_or_var_name+'_id')

        tiledb.from_pandas(
            uri=obs_or_var_uri,
            dataframe=obs_or_var_data,
            name=obs_or_var_name,
            sparse=True,
            allows_duplicates=False,
            offsets_filters=offsets_filters,
            attr_filters=attr_filters,
            dim_filters=dim_filters,
            capacity=100000,
            tile=extent,
            ctx=self.ctx
        )

        if self.verbose:
            print(util.format_elapsed(s, f"    FINISH WRITING {obs_or_var_uri}"))

        return obs_or_var_uri

    # ----------------------------------------------------------------
    def write_annotation_matrices(self, annotation_matrices, name: str, dim_name: str, dim_values):
        """
        Populates the obsm/ or varm/ subgroup for a SOMA object, then writes all the components
        arrays under that group.
        :param annotation_matrices: anndata.obsm or anndata.varm
        :param name: 'obsm' or 'varm'
        :param dim_name: 'obs_id' or 'var_id'
        :param dim_values: anndata.obs_names or anndata.var_names
        """
        assert name in ["obsm", "varm"]

        subgroup_uri = os.path.join(self.uri, name)
        tiledb.group_create(subgroup_uri, ctx=self.ctx)
        subgroup = tiledb.Group(subgroup_uri, mode="w", ctx=self.ctx)

        for mat_name in annotation_matrices.keys():
            mat = annotation_matrices[mat_name]
            component_array_uri = os.path.join(subgroup_uri, mat_name)
            if self.verbose:
                s = util.get_start_stamp()
                print(f"    START  WRITING {component_array_uri}")
                print(f"    Annotation matrix {name}/{mat_name} has shape {mat.shape}")

            # We do not have column names for anndata-provenance annotation matrices.
            # So, if say we're looking at anndata.obsm['X_pca'], we create column names
            # 'X_pca_1', 'X_pca_2', etc.
            (nrow, nattr) = mat.shape
            attr_names = [mat_name + '_' + str(j) for j in range(1, nattr+1)]

            # Ingest annotation matrices as 1D/multi-attribute sparse arrays
            self.__create_annot_matrix(component_array_uri, mat, mat_name, dim_name, attr_names)
            self.__ingest_annot_matrix(component_array_uri, mat, dim_values, attr_names)

            if self.verbose:
                print(util.format_elapsed(s, f"    FINISH WRITING {component_array_uri}"))

            subgroup.add(uri=component_array_uri, relative=False, name=mat_name)
        subgroup.close()

        return subgroup_uri

    # ----------------------------------------------------------------
    def write_annotation_pairwise_matrices(self, annotation_pairwise_matrices, name: str, dim_name: str, dim_values):
        """
        Populates the obsp/ or varp/ subgroup for a SOMA object, then writes all the components
        arrays under that group.
        :param annotation_matrices: anndata.obsp or anndata.varp
        :param name: 'obsp' or 'varp'
        :param dim_name: 'obs_id' or 'var_id'
        :param dim_values: anndata.obs_names or anndata.var_names
        """
        assert name in ["obsp", "varp"]

        subgroup_uri = os.path.join(self.uri, name)
        tiledb.group_create(subgroup_uri, ctx=self.ctx)
        subgroup = tiledb.Group(subgroup_uri, mode="w", ctx=self.ctx)

        dim_labels = [f"{dim_name}_{x}" for x in ["i", "j"]]

        for mat_name in annotation_pairwise_matrices.keys():
            mat = annotation_pairwise_matrices[mat_name]
            component_array_uri = os.path.join(subgroup_uri, mat_name)
            if self.verbose:
                s = util.get_start_stamp()
                print(f"    START  WRITING {component_array_uri}")
                print(f"    Annotation-pairwise matrix {name}/{mat_name} has shape {mat.shape}")

            # Ingest annotation pairwise matrices as 2D/single-attribute sparse arrays
            self.__create_coo_array(component_array_uri, dim_labels, "value", mat_dtype=mat.dtype)
            self.__ingest_coo_data_whole(component_array_uri, mat, dim_values, dim_values)

            if self.verbose:
                print(util.format_elapsed(s, f"    FINISH WRITING {component_array_uri}"))

            subgroup.add(uri=component_array_uri, relative=False, name=mat_name)
        subgroup.close()

        return subgroup_uri

    # ----------------------------------------------------------------
    def __create_annot_matrix(self, uri: str, mat, mat_name: str, dim_name: str, attr_names):
        """
        Create a TileDB 1D sparse array with string dimension and multiple attributes

        :param uri: URI of the array to be created
        :param mat: e.g. anndata.obsm['X_pca'] -- nominally a numpy.ndarray
        :param mat_name: e.g. 'X_pca' in anndata.obsm['X_pca']
        :param dim_name: name of the TileDB array dimension -- nominally 'obs_id' or 'var_id'
        :param attr_names: column names for the dataframe
        """
        assert isinstance(uri, str)
        assert isinstance(mat_name, str)
        assert isinstance(dim_name, str)

        # Nominally 'obs_id' or 'var_id'
        dom = tiledb.Domain(
            tiledb.Dim(name=dim_name, domain=(None, None), dtype="ascii", filters=[tiledb.ZstdFilter(level=22)]),
            ctx=self.ctx
        )

        # Verify:
        # anndata = ad.read_h5ad('anndata/pbmc3k_processed.h5ad')
        # anndata.obsm['X_pca'].dtype
        #dtype = 'float32'
        dtype = mat.dtype

        attrs = [
            tiledb.Attr(attr_name, dtype=dtype, filters=[tiledb.ZstdFilter()], ctx=self.ctx)
            for attr_name in attr_names
        ]

        sch = tiledb.ArraySchema(
            domain=dom,
            attrs=attrs,
            sparse=True,
            allows_duplicates=True,
            offsets_filters=[tiledb.DoubleDeltaFilter(), tiledb.BitWidthReductionFilter(), tiledb.ZstdFilter()],
            capacity=100000,
            cell_order='row-major',
            # As of TileDB core 2.8.2, we cannot consolidate string-indexed sparse arrays with
            # col-major tile order: so we write `X` with row-major tile order.
            tile_order='row-major',
            ctx=self.ctx
        )

        tiledb.Array.create(uri, sch, ctx=self.ctx)

    # ----------------------------------------------------------------
    def __ingest_annot_matrix(self, uri: str, mat, dim_values, col_names):
        """
        Convert ndarray/(csr|csc)matrix to a dataframe and ingest into TileDB.
        :param uri: TileDB URI of the array to be written.
        :param mat: Matrix-like object coercible to a pandas dataframe.
        :param dim_values: barcode/gene IDs from anndata.obs_names or anndata.var_names
        :param col_names: List of column names.
        """

        assert len(col_names) == mat.shape[1]

        df = pd.DataFrame(mat, columns = col_names)

        with tiledb.open(uri, mode="w", ctx=self.ctx) as A:
            A[dim_values] = df.to_dict(orient='list')

    # ----------------------------------------------------------------
    def __create_coo_array(self, uri: str, dim_labels, attr_name: str, mat_dtype):
        """
        Create a TileDB 2D sparse array with string dimensions and a single attribute.

        :param uri: URI of the array to be created
        :param mat: scipy.sparse.coo_matrix
        :param dim_labels: names of the TileDB array dimensions
        :param attr_name: name of the TileDB array-data attribute
        """
        assert isinstance(uri, str)
        assert len(dim_labels) == 2
        assert isinstance(attr_name, str)

        dom = tiledb.Domain(
            tiledb.Dim(name=dim_labels[0], domain=(None, None), dtype="ascii", filters=[tiledb.RleFilter()]),
            tiledb.Dim(name=dim_labels[1], domain=(None, None), dtype="ascii", filters=[tiledb.ZstdFilter(level=22)]),
            ctx=self.ctx
        )

        att = tiledb.Attr(attr_name, dtype=mat_dtype, filters=[tiledb.ZstdFilter()], ctx=self.ctx)
        sch = tiledb.ArraySchema(
            domain=dom,
            attrs=(att,),
            sparse=True,
            allows_duplicates=True,
            offsets_filters=[tiledb.DoubleDeltaFilter(), tiledb.BitWidthReductionFilter(), tiledb.ZstdFilter()],
            capacity=100000,
            cell_order='row-major',
            tile_order='col-major',
            ctx=self.ctx
        )
        tiledb.Array.create(uri, sch, ctx=self.ctx)

    # ----------------------------------------------------------------
    # Example: suppose this 4x3 is to be written in two chunks of two rows each
    # but written in sorted order.
    #
    # Original     Sorted     Permutation
    #  data       row names
    #
    #   X Y Z
    # C 0 1 2      A            1
    # A 4 0 5      B            2
    # B 7 0 0      C            0
    # D 0 8 9      D            3
    #
    # First chunk:
    # * Row indices 0,1 map to permutation indices 1,2
    # * i,i2 are 0,2
    # * chunk_coo is original matrix rows 1,2
    # * chunk_coo.row is [0,1]
    # * chunk_coo.row + i is [0,1]
    # * sorted_row_names: ['A', 'B']
    #
    # Second chunk:
    # * Row indices 2,3 map to permutation indices 0,3
    # * i,i2 are 2,4
    # * chunk_coo is original matrix rows 0,3
    # * chunk_coo.row is [0,1]
    # * chunk_coo.row + i is [2,3]
    # * sorted_row_names: ['C', 'D']
    #
    # See README-csr-ingest.md for important information of using this ingestor.
    # ----------------------------------------------------------------

    def __ingest_coo_data_rows_chunked(self, uri: str, mat: scipy.sparse._csr.csr_matrix, row_names, col_names):
        """
        Convert csr_matrix to coo_matrix chunkwise and ingest into TileDB.

        :param uri: TileDB URI of the array to be written.
        :param mat: csr_matrix.
        :param row_names: List of row names.
        :param col_names: List of column names.
        """

        assert len(row_names) == mat.shape[0]
        assert len(col_names) == mat.shape[1]

        # Sort the row names so we can write chunks indexed by sorted string keys.  This will lead
        # to efficient TileDB fragments in the sparse array indexed by these string keys.
        #
        # Key note: only the _obs labels_ are being sorted, and along with them come permutation
        # indices for accessing the CSR matrix via cursor-indirection -- e.g. csr[28] is accessed as
        # with csr[permuation[28]] -- the CSR matrix itself isn't sorted in bulk.
        sorted_row_names, permutation = util.get_sort_and_permutation(list(row_names))
        # Using numpy we can index this with a list of indices, which a plain Python list doesn't support.
        sorted_row_names = np.asarray(sorted_row_names)

        s = util.get_start_stamp()
        if self.verbose:
            print(f"      START  __ingest_coo_data_rows_chunked")

        with tiledb.open(uri, mode="w") as A:
            nrow = len(sorted_row_names)

            i = 0
            while i < nrow:
                # Find a number of CSR rows which will result in a desired nnz for the chunk.
                chunk_size = util.find_csr_chunk_size(mat, permutation, i, self.goal_chunk_nnz)
                i2 = i + chunk_size

                # Convert the chunk to a COO matrix.
                chunk_coo = mat[permutation[i:i2]].tocoo()

                s2 = util.get_start_stamp()

                # Write the chunk-COO to TileDB.
                d0 = sorted_row_names[chunk_coo.row + i]
                d1 = col_names[chunk_coo.col]

                if len(d0) == 0:
                    continue

                # Python ranges are (lo, hi) with lo inclusive and hi exclusive. But saying that
                # makes us look buggy if we say we're ingesting chunk 0:18 and then 18:32.
                # Instead, print inclusive lo..hi like 0..17 and 18..31.
                if self.verbose:
                    print("        START  chunk rows %d..%d of %d, obs_ids %s..%s, nnz=%d, %7.3f%%" %
                        (i, i2-1, nrow, d0[0], d0[-1], chunk_coo.nnz, 100*(i2-1)/nrow))

                # Write a TileDB fragment
                A[d0, d1] = chunk_coo.data

                if self.verbose:
                    print(util.format_elapsed(s2,"        FINISH chunk"))
                i = i2

        if self.verbose:
            print(util.format_elapsed(s,"      FINISH __ingest_coo_data_rows_chunked"))

    # ----------------------------------------------------------------
    def __ingest_coo_data_whole(self, uri: str, mat, row_names, col_names):
        """
        Convert ndarray/(csr|csc)matrix to coo_matrix and ingest into TileDB.

        :param uri: TileDB URI of the array to be written.
        :param mat: Matrix-like object coercible to a scipy coo_matrix.
        :param row_names: List of row names.
        :param col_names: List of column names.
        """

        assert len(row_names) == mat.shape[0]
        assert len(col_names) == mat.shape[1]

        mat_coo = scipy.sparse.coo_matrix(mat)
        d0 = row_names[mat_coo.row]
        d1 = col_names[mat_coo.col]

        with tiledb.open(uri, mode="w", ctx=self.ctx) as A:
            A[d0, d1] = mat_coo.data



    # ----------------------------------------------------------------
    def to_h5ad(self, h5ad_path: str):
        """
        Converts the soma group to anndata format and writes it to the specified .h5ad file.
        As of 2022-05-05 this is an incomplete prototype.
        """

        if self.verbose:
            s = util.get_start_stamp()
            print(f"START  SOMA.to_h5ad {self.uri} -> {h5ad_path}")

        anndata = self._to_anndata()

        if self.verbose:
            s2 = util.get_start_stamp()
            print(f"  START  write {h5ad_path}")
        anndata.write_h5ad(h5ad_path)
        if self.verbose:
            print(util.format_elapsed(s2, f"  FINISH write {h5ad_path}"))

        if self.verbose:
            print(util.format_elapsed(s, f"FINISH SOMA.to_h5ad {self.uri} -> {h5ad_path}"))

    # ----------------------------------------------------------------
    def to_anndata(self):
        """
        Converts the soma group to anndata. Choice of matrix formats is following
        what we often see in input .h5ad files:
        * X as scipy.sparse.csr_matrix
        * obs,var as pandas.dataframe
        * obsm,varm arrays as numpy.ndarray
        * obsp,varp arrays as scipy.sparse.csr_matrix
        As of 2022-05-05 this is an incomplete prototype.
        """

        if self.verbose:
            s = util.get_start_stamp()
            print(f"START  SOMA.to_anndata {self.uri}")

        retval = self._to_anndata()

        if self.verbose:
            print(util.format_elapsed(s, f"FINISH SOMA.to_anndata {self.uri}"))

        return retval

    # ----------------------------------------------------------------
    # Split from to_anndata solely to get the outermost verbosity prints done
    # exactly once, whether the user does soma.to_anndata or soma.to_h5ad.
    def _to_anndata(self):
        """
        Internal helper function for to_anndata; same arguments.
        """

        obs_df, obs_labels = self.outgest_obs_or_var('obs', 'obs_id')
        var_df, var_labels = self.outgest_obs_or_var('var', 'var_id')

        X_mat = self.outgest_X(obs_labels, var_labels)
        # TODO
        print("  X/RAW OUTGEST NOT IMPLMENTED YET")

        obsm = self.outgest_obsm_or_varm('obsm', 'obs_id')
        varm = self.outgest_obsm_or_varm('varm', 'var_id')

        # TODO
        print("  OBSP OUTGEST NOT WORKING YET")
        #obsp = self.outgest_obsp_or_varp('obsp')
        print("  VARP OUTGEST NOT WORKING YET")
        #varp = self.outgest_obsp_or_varp('varp')

        return ad.AnnData(
            X=X_mat, obs=obs_df, var=var_df, obsm=obsm, varm=varm,
        )

    # ----------------------------------------------------------------
    def outgest_obs_or_var(self, array_name: str, index_name: str):
        """
        Reads the TileDB obs or var array and returns a type of pandas dataframe
        and dimension values.
        :param array_name: 'obs' or 'var'
        :param index_name: 'obs_id' or 'var_id'
        """

        uri = f"{self.uri}/{array_name}"
        if self.verbose:
            s = util.get_start_stamp()
            print(f"  START  read {uri}")

        with tiledb.open(uri) as A:
            df = pd.DataFrame(A[:])
            labels = df[index_name] # strings, sorted
            df = df.set_index(index_name)
            retval = (df, labels)

        if self.verbose:
            print(util.format_elapsed(s, f"  FINISH read {uri}"))

        return retval

    # ----------------------------------------------------------------
    def outgest_X(self, obs_labels, var_labels):
        """
        Given a TileDB soma group, returns a scipy.sparse.csr_matrix with the X data.
        :param obs_labels: from the obs array. Note that TileDB will have sorted these.
        :param var_labels: from the var array. Note that TileDB will have sorted these.
        """

        uri = f"{self.uri}/X/data"
        if self.verbose:
            s = util.get_start_stamp()
            print(f"  START  read {uri}")

        # Since X is sparse, with two string dimensions, we get back a dict:
        # * 'obs_id' key is a sequence of dim0 coordinates for X data.
        # * 'var_id' key is a sequence of dim1 coordinates for X data.
        # * 'values' key is a sequence of X data values.
        with tiledb.open(uri) as X:
            X_data = X[:]

        # Now we need to convert from TileDB's string indices to CSR integer indices.
        # Make a dict from string dimension values to integer indices.
        #
        # Example: suppose the sparse matrix looks like:
        #
        #     S T U V
        #   A 4 . . 3
        #   B: 5 . 6 .
        #   C . 1 . 2
        #   D 8 7 . .
        #
        # The return value from the X[:] query is (obs_id,var_id,value) triples like
        #
        #   A,S,4 A,V,3 B,S,5 B,U,6 C,V,2 C,T,1 D,S,8 D,T,7
        #
        # whereas scipy csr is going to want
        #
        #   0,0,4 0,3,3 1,0,5 1,2,6 2,3,2 2,1,1 3,0,8 3,1,7
        #
        # In order to accomplish this, we need to map ['A','B','C','D'] to [0,1,2,3] via {'A':0,
        # 'B':1, 'C':2, 'D':3} and similarly for the other dimension.
        obs_labels_to_indices = dict(zip(obs_labels, [i for i,e in enumerate(obs_labels)]))
        var_labels_to_indices = dict(zip(var_labels, [i for i,e in enumerate(var_labels)]))

        # Apply the map.
        X_obs_indices = [obs_labels_to_indices[X_obs_label] for X_obs_label in X_data['obs_id']]
        X_var_indices = [var_labels_to_indices[X_var_label] for X_var_label in X_data['var_id']]

        retval = scipy.sparse.csr_matrix(
            (list(X_data['value']), (list(X_obs_indices), list(X_var_indices)))
        )

        if self.verbose:
            print(util.format_elapsed(s, f"  FINISH read {uri}"))

        return retval

    # ----------------------------------------------------------------
    def outgest_obsm_or_varm(self, group_name: str, index_name: str):
        """
        Reads the TileDB obsm or varm group and returns a dict from array name to array.
        These arrays are in numpy.ndarray format.
        :param soma_path: Path to the soma group.
        :param group_name: "obsm" or "varm".
        :param index_name: "obs_id" or "var_id".
        """
        group_uri = os.path.join(self.uri, group_name)

        grp = None
        try: # Not all groups have all four of obsm, obsp, varm, and varp.
            grp = tiledb.Group(group_uri, mode='r')
        except:
            pass
        if grp == None:
            if self.verbose:
                print(f"  {group_uri} not found")
            return {}

        if self.verbose:
            s = util.get_start_stamp()
            print(f"  START  read {group_uri}")

        matrices_in_group = {}
        for element in grp:
            with tiledb.open(element.uri) as A:
                with tiledb.open(element.uri) as A:
                    if self.verbose:
                        s2 = util.get_start_stamp()
                        print(f"    START  read {element.uri}")

                    df = pd.DataFrame(A[:])
                    df.set_index(index_name, inplace=True)
                    matrix_name = os.path.basename(element.uri) # e.g. 'X_pca'
                    matrices_in_group[matrix_name] = df.to_numpy()

                    if self.verbose:
                        print(util.format_elapsed(s2, f"    FINISH read {element.uri}"))

        grp.close()

        if self.verbose:
            print(util.format_elapsed(s, f"  FINISH read {group_uri}"))

        return matrices_in_group

    # ----------------------------------------------------------------
    def outgest_obsp_or_varp(self, group_name: str):
        """
        Reads the TileDB obsp or varp group and returns a dict from array name to array.
        These arrays are in scipy.csr format.
        :param soma_path: Path to the soma group.
        :param group_name: "obsp" or "varp".
        """
        group_uri = os.path.join(self.uri, group_name)

        grp = None
        try: # Not all groups have all four of obsm, obsp, varm, and varp.
            grp = tiledb.Group(group_uri, mode='r')
        except:
            pass
        if grp == None:
            if self.verbose:
                print(f"  {group_uri} not found")
            return {}

        if self.verbose:
            s = util.get_start_stamp()
            print(f"  START  read {group_uri}")

        matrices_in_group = {}
        for element in grp:
            with tiledb.open(element.uri) as A:
                with tiledb.open(element.uri) as A:
                    if self.verbose:
                        s2 = util.get_start_stamp()
                        print(f"    START  read {element.uri}")

                    df = pd.DataFrame(A[:])
                    matrix_name = os.path.basename(element.uri)
                    matrices_in_group[matrix_name] = scipy.sparse.coo_matrix(df).tocsr()
                    # TODO: not working yet:
                    # TypeError: no supported conversion for types: (dtype('O'),)

                    if self.verbose:
                        print(util.format_elapsed(s2, f"    FINISH read {element.uri}"))
        grp.close()

        if self.verbose:
            print(util.format_elapsed(s, f"  FINISH read {group_uri}"))

        return matrices_in_group
