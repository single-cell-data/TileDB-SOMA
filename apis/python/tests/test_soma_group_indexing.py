import anndata
import tiledb
import tiledbsc
import tiledbsc.io

import numpy as np

import pytest
import tempfile
import os
from pathlib import Path

HERE = Path(__file__).parent


@pytest.fixture
def h5ad_file(request):
    input_path = HERE.parent / "anndata/pbmc-small.h5ad"
    return input_path


def test_soma_group_indexing(h5ad_file):
    """
    Verify basic group-member access at the tiledbsc-py level.
    """

    # Set up anndata input path and tiledb-group output path
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name

    # Ingest
    soma = tiledbsc.SOMA(output_path, verbose=False)
    tiledbsc.io.from_h5ad(soma, h5ad_file)
    assert soma.exists()

    # Structure:
    #   X/data
    #   obs
    #   var
    #   obsm/X_pca
    #   obsm/X_tsne
    #   obsm/X_umap
    #   obsm/X_draw_graph_fr
    #   varm/PCs
    #   obsp/distances
    #   obsp/connectivities
    #   raw/X/data
    #   raw/var
    #   raw/varm/PCs

    assert set(soma._get_member_names()) == set(
        ["uns", "varm", "X", "raw", "obsp", "varp", "var", "obsm", "obs"]
    )
    assert set(soma.X._get_member_names()) == set(["data"])
    assert soma.X["data"].dim_names() == ["obs_id", "var_id"]
    assert soma.X["data"].shape() == (80, 20)

    assert soma.obs.exists()
    assert soma.obs.dim_names() == ["obs_id"]
    assert soma.obs.dim_name == "obs_id"
    assert soma.obs.keys() == [
        "orig.ident",
        "nCount_RNA",
        "nFeature_RNA",
        "RNA_snn_res.0.8",
        "letter.idents",
        "groups",
        "RNA_snn_res.1",
    ]
    assert soma.obs.shape() == (80, 7)
    assert set(soma.obs.ids()) == set(
        [
            "AAATTCGAATCACG",
            "AAGCAAGAGCTTAG",
            "AAGCGACTTTGACG",
            "AATGCGTGGACGGA",
            "AATGTTGACAGTCA",
            "ACAGGTACTGGTGT",
            "ACCAGTGAATACCG",
            "ACGTGATGCCATGA",
            "ACTCGCACGAAAGT",
            "AGAGATGATCTCGC",
            "AGATATACCCGTAA",
            "AGGTCATGAGTGTC",
            "AGTCAGACTGCACA",
            "AGTCTTACTTCGGA",
            "ATAAGTTGGTACGT",
            "ATACCACTCTAAGC",
            "ATAGGAGAAACAGA",
            "ATCATCTGACACCA",
            "ATGCCAGAACGACT",
            "ATTACCTGCCTTAT",
            "ATTCAGCTCATTGG",
            "ATTGCACTTGCTTT",
            "ATTGTAGATTCCCG",
            "CATATAGACTAAGC",
            "CATCAGGATGCACA",
            "CATCATACGGAGCA",
            "CATGAGACACGGGA",
            "CATGCGCTAGTCAC",
            "CATGGCCTGTGCAT",
            "CATTACACCAACTG",
            "CCATCCGATTCGCC",
            "CCCAACTGCAATCG",
            "CCTATAACGAGACG",
            "CGGCACGAACTCAG",
            "CGTAGCCTGTATGC",
            "CTAAACCTCTGACA",
            "CTAAACCTGTGCAT",
            "CTAACGGAACCGAT",
            "CTAGGTGATGGTTG",
            "CTGCCAACAGGAGC",
            "CTTCATGACCGAAT",
            "CTTGATTGATCTTC",
            "GAACCTGATGAACC",
            "GACATTCTCCACCT",
            "GACGCTCTCTCTCG",
            "GAGTTGTGGTAGCT",
            "GATAGAGAAGGGTG",
            "GATAGAGATCACGA",
            "GATATAACACGCAT",
            "GCACTAGACCTTTA",
            "GCAGCTCTGTTTCT",
            "GCGCACGACTTTAC",
            "GCGCATCTTGCTCC",
            "GCGTAAACACGGTT",
            "GCTCCATGAGAAGT",
            "GGAACACTTCAGAC",
            "GGCATATGCTTATC",
            "GGCATATGGGGAGT",
            "GGCCGATGTACTCT",
            "GGGTAACTCTAGTG",
            "GGTGGAGATTACTC",
            "GTAAGCACTCATTC",
            "GTCATACTTCGCCT",
            "GTTGACGATATCGG",
            "TACAATGATGCTAG",
            "TACATCACGCTAAC",
            "TACGCCACTCCGAA",
            "TACTCTGAATCGAC",
            "TAGGGACTGAACTC",
            "TCCACTCTGAGCTT",
            "TCTGATACACGTGT",
            "TGACTGGATTCTCA",
            "TGAGCTGAATGCTG",
            "TGGTATCTAAACAG",
            "TTACCATGAATCGC",
            "TTACGTACGTTCAG",
            "TTGAGGACTACGCA",
            "TTGCATTGAGCTAC",
            "TTGGTACTGAATCC",
            "TTTAGCTGTACTCT",
        ]
    )
    assert soma.obs.df().shape == (80, 7)
    assert soma.obs.df(["AAGCAAGAGCTTAG", "TTGGTACTGAATCC"]).shape == (2, 7)
    assert list(soma.obs.df().dtypes) == [
        np.dtype("int32"),
        np.dtype("float64"),
        np.dtype("int32"),
        np.dtype("int32"),
        np.dtype("int32"),
        np.dtype("O"),
        np.dtype("int32"),
    ]

    assert soma.var.exists()
    assert soma.var.dim_names() == ["var_id"]
    assert soma.obs.dim_name == "obs_id"
    assert soma.var.keys() == [
        "vst.mean",
        "vst.variance",
        "vst.variance.expected",
        "vst.variance.standardized",
        "vst.variable",
    ]
    assert soma.var.shape() == (20, 5)
    assert set(soma.var.ids()) == set(
        [
            "AKR1C3",
            "CA2",
            "CD1C",
            "GNLY",
            "HLA-DPB1",
            "HLA-DQA1",
            "IGLL5",
            "MYL9",
            "PARVB",
            "PF4",
            "PGRMC1",
            "PPBP",
            "RP11-290F20.3",
            "RUFY1",
            "S100A8",
            "S100A9",
            "SDPR",
            "TREML1",
            "TUBB1",
            "VDAC3",
        ]
    )
    assert soma.var.shape() == (20, 5)
    assert soma.var.df(["RUFY1", "AKR1C3"]).shape == (2, 5)
    assert list(soma.var.df().dtypes) == [
        np.dtype("float64"),
        np.dtype("float64"),
        np.dtype("float64"),
        np.dtype("float64"),
        np.dtype("int32"),
    ]

    assert soma.obsm.exists()
    assert set(soma.obsm._get_member_names()) == set(["X_pca", "X_tsne"])
    assert set(soma.obsm.keys()) == set(["X_pca", "X_tsne"])
    assert soma.obsm["X_pca"].exists()
    assert isinstance(soma.obsm["X_pca"], tiledbsc.AnnotationMatrix)
    assert soma.obsm["nonesuch"] is None
    assert soma.obsm["X_pca"].dim_names() == ["obs_id"]
    assert soma.obsm["X_pca"].shape() == (80, 19)
    assert soma.obsm["X_pca"].df().shape == (80, 19)
    assert list(soma.obsm["X_pca"].df().dtypes) == [
        np.dtype("float64"),
        np.dtype("float64"),
        np.dtype("float64"),
        np.dtype("float64"),
        np.dtype("float64"),
        np.dtype("float64"),
        np.dtype("float64"),
        np.dtype("float64"),
        np.dtype("float64"),
        np.dtype("float64"),
        np.dtype("float64"),
        np.dtype("float64"),
        np.dtype("float64"),
        np.dtype("float64"),
        np.dtype("float64"),
        np.dtype("float64"),
        np.dtype("float64"),
        np.dtype("float64"),
        np.dtype("float64"),
    ]

    assert set(soma.varm._get_member_names()) == set(["PCs"])
    assert soma.varm["PCs"].exists()
    assert isinstance(soma.varm["PCs"], tiledbsc.AnnotationMatrix)
    assert soma.varm["nonesuch"] is None
    assert soma.varm._get_member_names() == ["PCs"]
    assert soma.varm["PCs"].dim_names() == ["var_id"]
    assert soma.varm["PCs"].shape() == (20, 19)
    assert soma.varm["PCs"].df().shape == (20, 19)

    assert set(soma.obsp._get_member_names()) == set(["distances"])
    assert soma.obsp["distances"].exists()
    assert soma.obsp["distances"].dim_names() == ["obs_id_i", "obs_id_j"]
    assert soma.obsp["distances"].shape() == (80, 80)
    assert isinstance(soma.obsp["distances"], tiledbsc.AssayMatrix)

    assert soma.varp["nonesuch"] is None

    assert soma.uns.exists()
    assert set(soma.uns._get_member_names()) == set(["neighbors"])
    assert soma.uns["neighbors"].exists()
    assert soma.uns.exists()
    assert isinstance(soma.uns["neighbors"], tiledbsc.UnsGroup)
    assert set(soma.uns["neighbors"]._get_member_names()) == set(["params"])
    assert isinstance(soma.uns["neighbors"]["params"], tiledbsc.UnsGroup)
    assert set(soma.uns["neighbors"]["params"]._get_member_names()) == set(["method"])
    assert isinstance(soma.uns["neighbors"]["params"]["method"], tiledbsc.UnsArray)
    assert soma.uns["nonesuch"] is None

    assert list(soma.uns.keys()) == ["neighbors"]
    assert list(soma.uns["neighbors"].keys()) == ["params"]
    assert "neighbors" in soma.uns
    assert "nonesuch" not in soma.uns

    # We exercise these to make sure they're not throwing exceptions.
    for e in soma.obsm:
        foo = (e.name, e.df().shape, e.uri)
    for e in soma.varm:
        foo = (e.name, e.df().shape, e.uri)
    for e in soma.obsp:
        foo = (e.name, e.df().shape, e.uri)
    for e in soma.varp:
        foo = (e.name, e.df().shape, e.uri)


def test_not_exists():
    soma = tiledbsc.SOMA("/nonesuch/nowhere/never", verbose=False)
    assert not soma.exists()
    assert not soma.obs.exists()
