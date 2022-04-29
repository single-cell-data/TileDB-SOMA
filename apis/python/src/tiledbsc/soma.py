import os
from typing import Optional

import anndata as ad
import numpy   as np
import pandas  as pd
import pyarrow as pa
import scanpy
import scipy
import tiledb


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

        # If URI is "/something/test1" then:
        # * obs_uri  is "/something/test1/obs"
        # * var_uri  is "/something/test1/var"
        # * data_uri is "/something/test1/X"

        # If URI is "tiledb://namespace/s3://bucketname/something/test1" then:
        # * obs_uri  is "tiledb://namespace/s3://bucketname/something/test1/obs"
        # * var_uri  is "tiledb://namespace/s3://bucketname/something/test1/var"
        # * data_uri is "tiledb://namespace/s3://bucketname/something/test1/X"

    # ----------------------------------------------------------------
    def from_anndata(self, anndata: ad.AnnData):
        """
        Factory function to instantiate a SOMA object from an input anndata.AnnData object.
        """
        if self.verbose:
            print(f"START  SOMA.from_ann")

        anndata = self.decategoricalize(anndata)
        self.write_tiledb_group(anndata)

        if self.verbose:
            print("FINISH  SOMA.from_ann")

    # ----------------------------------------------------------------
    def from_h5ad(self, input_path: str):
        """
        Factory function to instantiate a SOMA object from an input .h5ad file.
        """
        if self.verbose:
            print(f"START  SOMA.from_h5ad {input_path} -> {self.uri}")

        anndata = self.read_h5ad(input_path)
        anndata = self.decategoricalize(anndata)
        self.write_tiledb_group(anndata)

        if self.verbose:
            print(f"FINISH SOMA.from_h5ad {input_path} -> {self.uri}")

    # ----------------------------------------------------------------
    def from_10x(self, input_path: str):
        """
        Factory function to instantiate a SOMA object from an input 10X file.
        """
        if self.verbose:
            print(f"START  SOMA.from_10x {input_path} -> {self.uri}")

        anndata = self.read_10x(input_path)
        anndata = self.decategoricalize(anndata)
        self.write_tiledb_group(anndata)

        if self.verbose:
            print(f"FINISH SOMA.from_10x {input_path} -> {self.uri}")

    # ----------------------------------------------------------------
    def read_h5ad(self, input_path: str):
        """
        File-ingestor for .h5ad files
        """
        if self.verbose:
            print(f"  START  READING {input_path}")
        anndata = ad.read_h5ad(input_path)
        anndata.var_names_make_unique()
        if self.verbose:
            print(f"  FINISH READING {input_path}")
        return anndata

    # ----------------------------------------------------------------
    def read_10x(self, input_path: str):
        """
        File-ingestor for 10X files
        """
        if self.verbose:
            print(f"  START  READING {input_path}")
        anndata = scanpy.read_10x_h5(input_path)
        anndata.var_names_make_unique()
        if self.verbose:
            print(f"  FINISH READING {input_path}")
        return anndata

    # ----------------------------------------------------------------
    def decategoricalize(self, anndata: ad.AnnData):
        """
        Performs a typecast for the categorical datatype that pandas can handle.
        We cast this to Python 'object' which in this case means 'string'.
        The uns and raw datasets are deferred.
        """

        if self.verbose:
            print(f"  START  DECATEGORICALIZING")

        # See also https://docs.scipy.org/doc/numpy-1.10.1/reference/arrays.dtypes.html
        uncat = lambda x: x.astype("O") if isinstance(x.dtype, pd.CategoricalDtype) else x

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
            raw=anndata.raw,  # expect Python 'None' type when there is no raw -- assignment OK
            obs=obs,
            var=var,
            obsm=anndata.obsm,
            obsp=anndata.obsp,
            varm=anndata.varm,
            varp=anndata.varp,
        )

        if self.verbose:
            print(f"  FINISH DECATEGORICALIZING")
        return anndata

    # ----------------------------------------------------------------
    def write_tiledb_group(self, anndata: ad.AnnData):
        """
        Top-level writer method for creating a TileDB group for a SOMA object.
        """
        if self.verbose:
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
            obsm_uri = self.write_annotation_matrices(anndata.obsm, "obsm", anndata.obs_names)
            base_group.add(uri=obsm_uri, relative=False, name="obsm")

        if len(anndata.varm.keys()) > 0:
            varm_uri = self.write_annotation_matrices(anndata.varm, "varm", anndata.var_names)
            base_group.add(uri=varm_uri, relative=False, name="varm")

        if len(anndata.obsp.keys()) > 0:
            obsp_uri = self.write_annotation_pairwise_matrices(anndata.obsp, "obsp", "obs", anndata.obs_names)
            base_group.add(uri=obsp_uri, relative=False, name="obsp")

        if len(anndata.varp.keys()) > 0:
            varp_uri = self.write_annotation_pairwise_matrices(anndata.varp, "varp", "var", anndata.var_names)
            base_group.add(uri=varp_uri, relative=False, name="varp")

        # ----------------------------------------------------------------
        if self.verbose:
            print(f"  FINISH WRITING {self.uri}")

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
        * x is anndata.X or raw
        * group_uri is the URI of the parent group, e.g. 'foo/X'
        * arrayname is the name of the array within the parent group, e.g. 'data' for 'foo/X/data'
        * obs_names and var_names are the names for the axes
        """
        X_array_uri = os.path.join(group_uri, arrayname)
        if self.verbose:
            print(f"    START  WRITING {X_array_uri}")

        # Here we do not use tiledb.from_numpy, so that we can have more control over the schema.
        obs_dim, var_dim = np.meshgrid(obs_names, var_names)

        self.__create_2d_coo_array(uri=X_array_uri, dim_names=["obs_id", "var_id"], attr_name="value")
        self.__ingest_2d_coo_data(X_array_uri, x, obs_names, var_names)

        if self.verbose:
            print(f"    FINISH WRITING {X_array_uri}")
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
            print(f"    FINISH WRITING {obs_or_var_uri}")

        return obs_or_var_uri

    # ----------------------------------------------------------------
    def write_annotation_matrices(self, annotation_matrices, name: str, dim_names):
        """
        Populates the obsm/ or varm/ subgroup for a SOMA object, then writes all the components
        arrays under that group.
        * annotation_matrices: anndata.obsm or anndata.varm
        * name: 'obsm' or 'varm'
        * dim_names: anndata.obs_names or anndata.var_names
        """
        assert name in ["obsm", "varm"]

        subgroup_uri = os.path.join(self.uri, name)
        tiledb.group_create(subgroup_uri, ctx=self.ctx)
        subgroup = tiledb.Group(subgroup_uri, mode="w", ctx=self.ctx)

        # dim_labels = [f"{name}_id_{x}" for x in ["i", "j"]]

        for mat_name in annotation_matrices.keys():
            mat = annotation_matrices[mat_name]
            print(f"    Annotation matrix {name}/{mat_name} has shape {mat.shape}")
            component_array_uri = os.path.join(subgroup_uri, mat_name)
            if self.verbose:
                print(f"    START  WRITING {component_array_uri}")

            # Ingest annotation matrices as 1D sparse arrays
            # **TODO**

            # Check for conversion from pandas if necessary.  For the pbmc3k_processed reference
            # dataset, obsm and varm matrices are numpy.ndarray while obsp matrices are
            # scipy.sparse.csr.csr_matrix. For ongoing work we will likely need more checks
            # here. See also desc-ann.py in this directory which helps reveal the datatypes
            # contained within a given HDF5 file.
            if isinstance(mat, scipy.sparse.csr_matrix):
                mat = mat.toarray()

            tiledb.from_numpy(
                uri=component_array_uri,
                array=mat,
                ctx=self.ctx
            )

            if self.verbose:
                print(f"    FINISH WRITING {component_array_uri}")

            subgroup.add(uri=component_array_uri, relative=False, name=mat_name)
        subgroup.close()

        return subgroup_uri

    # ----------------------------------------------------------------
    def write_annotation_pairwise_matrices(self, annotation_pairwise_matrices, name: str, dim_names_name: str, dim_names):
        """
        Populates the obsp/ or varp/ subgroup for a SOMA object, then writes all the components
        arrays under that group.
        * annotation_matrices: anndata.obsp or anndata.varp
        * name: 'obsp' or 'varp'
        * dim_names_name: 'obs' or 'var'
        * dim_names: anndata.obs_names or anndata.var_names
        """
        assert name in ["obsp", "varp"]

        subgroup_uri = os.path.join(self.uri, name)
        tiledb.group_create(subgroup_uri, ctx=self.ctx)
        subgroup = tiledb.Group(subgroup_uri, mode="w", ctx=self.ctx)

        dim_labels = [f"{dim_names_name}_id_{x}" for x in ["i", "j"]]

        for mat_name in annotation_pairwise_matrices.keys():
            mat = annotation_pairwise_matrices[mat_name]
            print(f"    Annotation matrix {name}/{mat_name} has shape {mat.shape}")
            component_array_uri = os.path.join(subgroup_uri, mat_name)
            if self.verbose:
                print(f"    START  WRITING {component_array_uri}")

            # Ingest annotation pairwise matrices as 2D sparse arrays
            self.__create_2d_coo_array(component_array_uri, dim_labels, "value")
            self.__ingest_2d_coo_data(component_array_uri, mat, dim_names, dim_names)

            if self.verbose:
                print(f"    FINISH WRITING {component_array_uri}")

            subgroup.add(uri=component_array_uri, relative=False, name=mat_name)
        subgroup.close()

        return subgroup_uri

    # ----------------------------------------------------------------
    def __create_2d_coo_array(self, uri: str, dim_names, attr_name: str):
        """
        Create a TileDB 2D sparse array with string dimensions and a single attribute.

        :param uri: URI of the array to be created
        :param x: scipy.sparse.coo_matrix
        :param dim_names: names of the TileDB array dimensions
        :param attr_name: name of the TileDB array attribute
        """
        assert isinstance(uri, str)
        assert len(dim_names) == 2
        assert isinstance(attr_name, str)

        dom = tiledb.Domain(
            tiledb.Dim(name=dim_names[0], domain=(None, None), dtype="ascii", filters=[tiledb.RleFilter()]),
            tiledb.Dim(name=dim_names[1], domain=(None, None), dtype="ascii", filters=[tiledb.ZstdFilter(level=22)]),
            ctx=self.ctx
        )

        # Verify:
        # anndata = ad.read_h5ad('anndata/pbmc3k_processed.h5ad')
        # anndata.X.dtype
        dtype = 'float32'

        att = tiledb.Attr(attr_name, dtype=dtype, filters=[tiledb.ZstdFilter()], ctx=self.ctx)
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
    def __ingest_2d_coo_data(self, uri: str, x, row_names, col_names):
        """
        Convert ndarray/(csr|csc)matrix to coo_matrix and ingest into TileDB.

        :param uri: TileDB URI of the array to be written.
        :param x: Matrix-like object coercible to a scipy coo_matrix.
        :param row_names: List of row names.
        :param col_names: List of column names.
        """

        assert len(row_names) == x.shape[0]
        assert len(col_names) == x.shape[1]

        x_coo = scipy.sparse.coo_matrix(x)
        d0 = row_names[x_coo.row]
        d1 = col_names[x_coo.col]

        with tiledb.open(uri, mode="w", ctx=self.ctx) as A:
            A[d0, d1] = x_coo.data
