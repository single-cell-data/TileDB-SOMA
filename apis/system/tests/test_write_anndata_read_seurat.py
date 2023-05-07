import anndata as ad
import pytest

from tiledbsoma import io

from .common import TestWritePythonReadR, embed_python_list_into_R_code


class TestAnndataSOMASeurat(TestWritePythonReadR):
    """
    This test will perform an AnnData->SOMA->Seurat conversion and verify that the final artifact (Seurat)
    has the same properties of the initial one (AnnData)
    """

    @pytest.fixture(scope="class")
    def h5ad(self):
        """
        Fixture that will load an h5ad, convert it to SOMA, deploy it to `self.uri` and return a parsed version of it.
        """
        h5ad_path = "apis/python/testdata/pbmc-small.h5ad"
        io.from_h5ad(self.uri, input_path=h5ad_path, measurement_name="RNA")
        return ad.read_h5ad(h5ad_path)

    # Prepares an R script with the dependencies and loads the object in seuratObject
    # Future assertions can be done by appending to this script
    def base_R_script(self):
        return f"""
        library("tiledbsoma")
        library("SeuratObject")
        soma <- SOMAOpen("{self.uri}")
        ms <- soma$ms$get("RNA")
        query <- SOMAExperimentAxisQuery$new(
            experiment = soma,
            measurement_name = "RNA"
        )
        seuratObject <- query$to_seurat(c(data = "data"), var_index="var_id", obs_index="obs_id")
        """

    # tests
    def test_seurat_dim_matches(self, h5ad):
        n_var = len(h5ad.var)
        n_obs = len(h5ad.obs)
        self.r_assert(
            f"""
            stopifnot(ncol(seuratObject) == {n_obs})
            stopifnot(nrow(seuratObject) == {n_var})
        """
        )

    def test_seurat_var_names_match(self, h5ad):
        var_names = h5ad.var.index.tolist()
        var_names_list = embed_python_list_into_R_code(var_names)
        self.r_assert(
            f"""
            stopifnot(all.equal(rownames(x = seuratObject), c({var_names_list})))
        """
        )

    def test_seurat_obs_names_match(self, h5ad):
        obs_names = h5ad.obs.index.tolist()
        obs_names_list = embed_python_list_into_R_code(obs_names)
        self.r_assert(
            f"""
            stopifnot(all.equal(colnames(x = seuratObject), c({obs_names_list})))
        """
        )
