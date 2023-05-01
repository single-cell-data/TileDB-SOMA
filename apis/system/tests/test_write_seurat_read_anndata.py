import pytest

import tiledbsoma as soma
from tiledbsoma import io

from .common import (
    TestWritePythonReadR,
    embed_python_list_into_R_code,
)


class TestSeuratSOMAAnndata(TestWritePythonReadR):
    """
    This test will perform a Seurat->SOMA->Anndata conversion and verify that the final artifact (Anndata)
    has the same properties of the initial one (Seurat)
    """

    @pytest.fixture(scope="class")
    def anndata_from_seurat(self):
        base_script = f"""
        library("tiledbsoma")
        library("SeuratObject")

        data("pbmc_small")
        write_soma(pbmc_small, "{self.uri}")
        """
        self.execute_R_script(base_script)
        pbmc = soma.open(self.uri)
        ad = io.to_anndata(pbmc, measurement_name="RNA")
        return ad

    def base_R_script(self):
        # Since we can't read Seurat in python, we need to do assertions in R, hence this script
        return """
        library("SeuratObject")
        data("pbmc_small")
        seuratObject <- pbmc_small
        """

    # tests
    def test_anndata_dim_matches(self, anndata_from_seurat):
        n_var = len(anndata_from_seurat.var)
        n_obs = len(anndata_from_seurat.obs)
        self.r_assert(
            f"""
            stopifnot(ncol(seuratObject) == {n_obs})
            stopifnot(nrow(seuratObject) == {n_var})
            """
        )

    def test_anndata_var_names_match(self, anndata_from_seurat):
        var_names = anndata_from_seurat.var.index.tolist()
        var_names_list = embed_python_list_into_R_code(var_names)
        self.r_assert(
            f"""
            stopifnot(all.equal(rownames(x = seuratObject), c({var_names_list})))
            """
        )

    def test_seurat_obs_names_match(self, anndata_from_seurat):
        obs_names = anndata_from_seurat.obs.index.tolist()
        obs_names_list = embed_python_list_into_R_code(obs_names)
        self.r_assert(
            f"""
            stopifnot(all.equal(colnames(x = seuratObject), c({obs_names_list})))
        """
        )
