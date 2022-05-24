import os
from typing import Optional, Union, Dict

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow as pa
import scanpy
import scipy
import tiledb

from tiledbsc import util
from tiledbsc import util_ann

from .soma_options import SOMAOptions
from .tiledb_group import TileDBGroup
from .assay_matrix_group import AssayMatrixGroup
from .annotation_dataframe import AnnotationDataFrame
from .annotation_matrix_group import AnnotationMatrixGroup
from .annotation_pairwise_matrix_group import AnnotationPairwiseMatrixGroup
from .raw_group import RawGroup
from .uns_group import UnsGroup


class SOMA(TileDBGroup):
    """Single-cell group
    Class for representing a group of TileDB groups/arrays that constitute an SOMA ('stack of matrices, annotated')
    which includes:
    * `X` (`AssayMatrixGroup`): a group of one or more labeled 2D sparse arrays
      that share the same dimensions.
    * `obs` (`AnnotationDataframe`): 1D labeled array with column labels for
      `X`
    * `var` ([`AnnotationDataframe`]): 1D labeled array with row labels for `X`

    See also desc-ann.py in this directory for helpful information to
    reveal the diversity/variety of HDF5 files we process.
    """

    X: AssayMatrixGroup
    obs: AnnotationDataFrame
    var: AnnotationDataFrame
    obsm: AnnotationMatrixGroup
    varm: AnnotationMatrixGroup
    obsp: AnnotationPairwiseMatrixGroup
    varp: AnnotationPairwiseMatrixGroup
    raw: RawGroup
    uns: UnsGroup

    # ----------------------------------------------------------------
    def __init__(
        self,
        uri: str,
        name="soma",
        soma_options: Optional[SOMAOptions] = None,
        verbose: Optional[bool] = True,
        config: Optional[tiledb.Config] = None,
        ctx: Optional[tiledb.Ctx] = None,
        parent: Optional[TileDBGroup] = None,  # E.g. a SOMA collection
    ):
        """
        @description Create a new SOMA object. The existing array group is
          opened at the specified array `uri` if one is present, otherwise a new
          array group is created.
        @param uri URI of the TileDB group
        @param verbose Print status messages
        """

        if ctx is None and config is not None:
            ctx = tiledb.Ctx(config)
        if soma_options is None:
            soma_options = SOMAOptions()  # Use default values from the constructor
        super().__init__(
            uri=uri,
            name=name,
            parent=parent,
            verbose=verbose,
            soma_options=soma_options,
            ctx=ctx,
        )

        X_uri = os.path.join(self.uri, "X")
        obs_uri = os.path.join(self.uri, "obs")
        var_uri = os.path.join(self.uri, "var")
        obsm_uri = os.path.join(self.uri, "obsm")
        varm_uri = os.path.join(self.uri, "varm")
        obsp_uri = os.path.join(self.uri, "obsp")
        varp_uri = os.path.join(self.uri, "varp")
        raw_uri = os.path.join(self.uri, "raw")
        uns_uri = os.path.join(self.uri, "uns")

        self.X = AssayMatrixGroup(
            uri=X_uri,
            name="X",
            row_dim_name="obs_id",
            col_dim_name="var_id",
            parent=self,
        )
        self.obs = AnnotationDataFrame(uri=obs_uri, name="obs", parent=self)
        self.var = AnnotationDataFrame(uri=var_uri, name="var", parent=self)
        self.obsm = AnnotationMatrixGroup(uri=obsm_uri, name="obsm", parent=self)
        self.varm = AnnotationMatrixGroup(uri=varm_uri, name="varm", parent=self)
        self.obsp = AnnotationPairwiseMatrixGroup(
            uri=obsp_uri, name="obsp", parent=self
        )
        self.varp = AnnotationPairwiseMatrixGroup(
            uri=varp_uri, name="varp", parent=self
        )
        self.raw = RawGroup(uri=raw_uri, name="raw", parent=self)
        self.uns = UnsGroup(uri=uns_uri, name="uns", parent=self)

        # If URI is "/something/test1" then:
        # * obs_uri  is "/something/test1/obs"
        # * var_uri  is "/something/test1/var"
        # * data_uri is "/something/test1/X"

        # If URI is "tiledb://namespace/s3://bucketname/something/test1" then:
        # * obs_uri  is "tiledb://namespace/s3://bucketname/something/test1/obs"
        # * var_uri  is "tiledb://namespace/s3://bucketname/something/test1/var"
        # * data_uri is "tiledb://namespace/s3://bucketname/something/test1/X"

    # ----------------------------------------------------------------
    def __str__(self):
        """
        Implements 'print(soma)'.
        """
        return f"name={self.name},uri={self.uri}"

    # ----------------------------------------------------------------
    def from_h5ad(self, input_path: str):
        """
        Reads an .h5ad file and writes to a TileDB group structure.
        """
        if self._verbose:
            s = util.get_start_stamp()
            print(f"START  SOMA.from_h5ad {input_path} -> {self.uri}")

        anndata = self.read_h5ad(input_path)
        self.from_anndata(anndata)

        if self._verbose:
            print(
                util.format_elapsed(
                    s, f"FINISH SOMA.from_h5ad {input_path} -> {self.uri}"
                )
            )

    # ----------------------------------------------------------------
    def from_10x(self, input_path: str):
        """
        Reads a 10X file and writes to a TileDB group structure.
        """
        if self._verbose:
            s = util.get_start_stamp()
            print(f"START  SOMA.from_10x {input_path} -> {self.uri}")

        anndata = self.read_10x(input_path)

        self.from_anndata(anndata)

        if self._verbose:
            print(
                util.format_elapsed(
                    s, f"FINISH SOMA.from_10x {input_path} -> {self.uri}"
                )
            )

    # ----------------------------------------------------------------
    def read_h5ad(self, input_path: str):
        """
        File-ingestor for .h5ad files
        """
        if self._verbose:
            s = util.get_start_stamp()
            print(f"{self._indent}START  READING {input_path}")

        anndata = ad.read_h5ad(input_path)

        if self._verbose:
            print(util.format_elapsed(s, f"{self._indent}FINISH READING {input_path}"))
        return anndata

    # ----------------------------------------------------------------
    def read_10x(self, input_path: str):
        """
        File-ingestor for 10X files
        """
        if self._verbose:
            s = util.get_start_stamp()
            print(f"{self._indent}START  READING {input_path}")

        anndata = scanpy.read_10x_h5(input_path)

        if self._verbose:
            print(util.format_elapsed(s, f"{self._indent}FINISH READING {input_path}"))
        return anndata

    # ================================================================
    # WRITE PATH INTO TILEDB

    def from_anndata(self, anndata: ad.AnnData):
        """
        Top-level writer method for creating a TileDB group for a SOMA object.
        """

        # Without _at least_ an index, there is nothing to indicate the dimension indices.
        if anndata.obs.index.empty or anndata.var.index.empty:
            raise NotImplementedError("Empty AnnData.obs or AnnData.var unsupported.")

        if self._verbose:
            s = util.get_start_stamp()
            print(f"{self._indent}START  DECATEGORICALIZING")
        anndata.obs_names_make_unique()
        anndata.var_names_make_unique()
        anndata = util_ann._decategoricalize(anndata)
        if self._verbose:
            print(util.format_elapsed(s, f"{self._indent}FINISH DECATEGORICALIZING"))

        if self._verbose:
            s = util.get_start_stamp()
            print(f"{self._indent}START  WRITING {self.uri}")

        # Must be done first, to create the parent directory
        with self._open("w") as G:

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            self.X.from_matrix(anndata.X, anndata.obs.index, anndata.var.index)
            self._add_object(G, self.X)

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            self.obs.from_dataframe(dataframe=anndata.obs, extent=256)
            self._add_object(G, self.obs)

            self.var.from_dataframe(dataframe=anndata.var, extent=2048)
            self._add_object(G, self.var)

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            self.obsm.from_anndata(anndata.obsm, anndata.obs_names)
            self._add_object(G, self.obsm)

            self.varm.from_anndata(anndata.varm, anndata.var_names)
            self._add_object(G, self.varm)

            self.obsp.from_anndata(anndata.obsp, anndata.obs_names)
            self._add_object(G, self.obsp)

            self.varp.from_anndata(anndata.varp, anndata.var_names)
            self._add_object(G, self.varp)

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            if anndata.raw != None:
                self.raw.from_anndata(anndata)
                self._add_object(G, self.raw)

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            if anndata.uns != None:
                self.uns.from_anndata_uns(anndata.uns)
                self._add_object(G, self.uns)

        if self._verbose:
            print(util.format_elapsed(s, f"{self._indent}FINISH WRITING {self.uri}"))

    # ================================================================
    # READ PATH OUT OF TILEDB

    # ----------------------------------------------------------------
    def to_h5ad(self, h5ad_path: str):
        """
        Converts the soma group to anndata format and writes it to the specified .h5ad file.
        As of 2022-05-05 this is an incomplete prototype.
        """

        if self._verbose:
            s = util.get_start_stamp()
            print(f"START  SOMA.to_h5ad {self.uri} -> {h5ad_path}")

        anndata = self._to_anndata()

        if self._verbose:
            s2 = util.get_start_stamp()
            print(f"{self._indent}START  write {h5ad_path}")
        anndata.write_h5ad(h5ad_path)
        if self._verbose:
            print(util.format_elapsed(s2, f"{self._indent}FINISH write {h5ad_path}"))

        if self._verbose:
            print(
                util.format_elapsed(s, f"FINISH SOMA.to_h5ad {self.uri} -> {h5ad_path}")
            )

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

        if self._verbose:
            s = util.get_start_stamp()
            print(f"START  SOMA.to_anndata {self.uri}")

        retval = self._to_anndata()

        if self._verbose:
            print(util.format_elapsed(s, f"FINISH SOMA.to_anndata {self.uri}"))

        return retval

    # ----------------------------------------------------------------
    # Split from to_anndata solely to get the outermost verbosity prints done
    # exactly once, whether the user does soma.to_anndata or soma.to_h5ad.
    def _to_anndata(self) -> ad.AnnData:
        """
        Internal helper function for to_anndata; same arguments.
        Note: obsp/varp outgest is omitted at present, pending some issues.
        """

        obs_df = self.obs.to_dataframe()
        var_df = self.var.to_dataframe()

        X_mat = self.X.data.to_csr_matrix(obs_df.index, var_df.index)

        obsm = self.obsm.to_dict_of_csr()
        varm = self.varm.to_dict_of_csr()

        # TODO
        print("  OBSP OUTGEST NOT WORKING YET")
        # obsp = self.obsp.to_dict_of_csr()
        print("  VARP OUTGEST NOT WORKING YET")
        # varp = self.varp.to_dict_of_csr()

        (raw_X, raw_var_df, raw_varm) = self.raw.to_anndata_raw(obs_df.index)

        anndata = ad.AnnData(X=X_mat, obs=obs_df, var=var_df, obsm=obsm, varm=varm)

        raw = ad.Raw(
            anndata,
            X=raw_X,
            var=raw_var_df,
            varm=raw_varm,
        )

        uns = self.uns.to_dict_of_matrices()

        anndata = ad.AnnData(
            X=anndata.X,
            dtype=None
            if anndata.X is None
            else anndata.X.dtype,  # some datasets have no X
            obs=anndata.obs,
            var=anndata.var,
            obsm=anndata.obsm,
            obsp=anndata.obsp,
            varm=anndata.varm,
            varp=anndata.varp,
            raw=raw,
            uns=uns,
        )

        return anndata

    # ----------------------------------------------------------------
    def to_anndata_from_raw(self):
        """
        Extract only the raw parts as a new AnnData object.
        """

        obs_df = self.obs.to_dataframe()
        var_df = self.raw.var.to_dataframe()
        X_mat = self.raw.X.data.to_csr_matrix(obs.index, var.index)

        return ad.AnnData(
            X=X_mat,
            obs=obs_df,
            var=var_df,
        )
