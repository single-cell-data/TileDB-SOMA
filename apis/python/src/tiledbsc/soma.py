import os
from typing import Optional

import anndata as ad
import numpy   as np
import pandas  as pd
import pyarrow as pa
import scanpy
import scipy
import tiledb
import tiledbsc.util     as util

def uncat(x):
    if isinstance(x.dtype, pd.CategoricalDtype):
        return x.astype('O')
    elif x.dtype == 'bool':
        return x.astype('uint8')
    elif x.dtype == np.float16:
        return x.astype(np.float16)
    else:
        return x

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
        anndata.var_names_make_unique()
        if self.verbose:
            print(util.format_elapsed(s, f"  FINISH READING {input_path}"))
        return anndata

    # ----------------------------------------------------------------
    def decategoricalize(self, anndata: ad.AnnData):
        """
        Performs a typecast for the categorical datatype that pandas can handle.
        We cast this to Python 'object' which in this case means 'string'.
        The uns and raw datasets are deferred.
        """

        if self.verbose:
            s = util.get_start_stamp()
            print(f"  START  DECATEGORICALIZING")

        # See also https://docs.scipy.org/doc/numpy-1.10.1/reference/arrays.dtypes.html
        #uncat = lambda x: x.astype("O") if isinstance(x.dtype, pd.CategoricalDtype) else x

        obs = pd.DataFrame.from_dict({k: uncat(v) for k, v in anndata.obs.items()})
        var = pd.DataFrame.from_dict({k: uncat(v) for k, v in anndata.var.items()})

        for key in anndata.obsm.keys():
            anndata.obsm[key] = uncat(anndata.obsm[key])
        for key in anndata.varm.keys():
            anndata.varm[key] = uncat(anndata.varm[key])
        for key in anndata.obsp.keys():
            anndata.obsp[key] = uncat(anndata.obsp[key])
        for key in anndata.varp.keys():
            anndata.varp[key] = uncat(anndata.varp[key])

        anndata = ad.AnnData(
            X=anndata.X,
            dtype=anndata.X.dtype,
            raw=anndata.raw,  # expect Python 'None' type when there is no raw -- assignment OK
            obs=obs,
            var=var,
            obsm=anndata.obsm,
            obsp=anndata.obsp,
            varm=anndata.varm,
            varp=anndata.varp,
        )

        if self.verbose:
            print(util.format_elapsed(s, f"  FINISH DECATEGORICALIZING"))
        return anndata

    # ----------------------------------------------------------------
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
                print("        START  chunk rows %d..%d of %d, obs_ids %s..%s, nnz=%d, %7.3f%%" %
                    (i, i2-1, nrow, d0[0], d0[-1], chunk_coo.nnz, 100*(i2-1)/nrow))

                # Write a TileDB fragment
                A[d0, d1] = chunk_coo.data

                print(util.format_elapsed(s2,"        FINISH chunk"))
                i = i2

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
