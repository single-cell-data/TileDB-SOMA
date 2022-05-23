import anndata
import tiledb
import tiledbsc

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
    soma.from_h5ad(h5ad_file)
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
    assert soma.X.data.dim_names() == ["obs_id", "var_id"]

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
    assert set(soma.obs.ids()) == set(
        [
            b"AAATTCGAATCACG",
            b"AAGCAAGAGCTTAG",
            b"AAGCGACTTTGACG",
            b"AATGCGTGGACGGA",
            b"AATGTTGACAGTCA",
            b"ACAGGTACTGGTGT",
            b"ACCAGTGAATACCG",
            b"ACGTGATGCCATGA",
            b"ACTCGCACGAAAGT",
            b"AGAGATGATCTCGC",
            b"AGATATACCCGTAA",
            b"AGGTCATGAGTGTC",
            b"AGTCAGACTGCACA",
            b"AGTCTTACTTCGGA",
            b"ATAAGTTGGTACGT",
            b"ATACCACTCTAAGC",
            b"ATAGGAGAAACAGA",
            b"ATCATCTGACACCA",
            b"ATGCCAGAACGACT",
            b"ATTACCTGCCTTAT",
            b"ATTCAGCTCATTGG",
            b"ATTGCACTTGCTTT",
            b"ATTGTAGATTCCCG",
            b"CATATAGACTAAGC",
            b"CATCAGGATGCACA",
            b"CATCATACGGAGCA",
            b"CATGAGACACGGGA",
            b"CATGCGCTAGTCAC",
            b"CATGGCCTGTGCAT",
            b"CATTACACCAACTG",
            b"CCATCCGATTCGCC",
            b"CCCAACTGCAATCG",
            b"CCTATAACGAGACG",
            b"CGGCACGAACTCAG",
            b"CGTAGCCTGTATGC",
            b"CTAAACCTCTGACA",
            b"CTAAACCTGTGCAT",
            b"CTAACGGAACCGAT",
            b"CTAGGTGATGGTTG",
            b"CTGCCAACAGGAGC",
            b"CTTCATGACCGAAT",
            b"CTTGATTGATCTTC",
            b"GAACCTGATGAACC",
            b"GACATTCTCCACCT",
            b"GACGCTCTCTCTCG",
            b"GAGTTGTGGTAGCT",
            b"GATAGAGAAGGGTG",
            b"GATAGAGATCACGA",
            b"GATATAACACGCAT",
            b"GCACTAGACCTTTA",
            b"GCAGCTCTGTTTCT",
            b"GCGCACGACTTTAC",
            b"GCGCATCTTGCTCC",
            b"GCGTAAACACGGTT",
            b"GCTCCATGAGAAGT",
            b"GGAACACTTCAGAC",
            b"GGCATATGCTTATC",
            b"GGCATATGGGGAGT",
            b"GGCCGATGTACTCT",
            b"GGGTAACTCTAGTG",
            b"GGTGGAGATTACTC",
            b"GTAAGCACTCATTC",
            b"GTCATACTTCGCCT",
            b"GTTGACGATATCGG",
            b"TACAATGATGCTAG",
            b"TACATCACGCTAAC",
            b"TACGCCACTCCGAA",
            b"TACTCTGAATCGAC",
            b"TAGGGACTGAACTC",
            b"TCCACTCTGAGCTT",
            b"TCTGATACACGTGT",
            b"TGACTGGATTCTCA",
            b"TGAGCTGAATGCTG",
            b"TGGTATCTAAACAG",
            b"TTACCATGAATCGC",
            b"TTACGTACGTTCAG",
            b"TTGAGGACTACGCA",
            b"TTGCATTGAGCTAC",
            b"TTGGTACTGAATCC",
            b"TTTAGCTGTACTCT",
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
    assert set(soma.var.ids()) == set(
        [
            b"AKR1C3",
            b"CA2",
            b"CD1C",
            b"GNLY",
            b"HLA-DPB1",
            b"HLA-DQA1",
            b"IGLL5",
            b"MYL9",
            b"PARVB",
            b"PF4",
            b"PGRMC1",
            b"PPBP",
            b"RP11-290F20.3",
            b"RUFY1",
            b"S100A8",
            b"S100A9",
            b"SDPR",
            b"TREML1",
            b"TUBB1",
            b"VDAC3",
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
    assert soma.obsm["X_pca"].df().shape == (80, 20)
    assert list(soma.obsm["X_pca"].df().dtypes) == [
        np.dtype("O"),
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

    assert set(soma.obsp._get_member_names()) == set(["distances"])
    assert soma.obsp["distances"].exists()
    assert soma.obsp["distances"].dim_names() == ["obs_id_i", "obs_id_j"]
    assert isinstance(soma.obsp["distances"], tiledbsc.AnnotationPairwiseMatrix)

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

    # These print to stdout -- here we exercise just them to make sure they're not throwing
    # exceptions.
    for e in soma.obsm:
        print(e.name, e.shape(), e.uri)
    for e in soma.varm:
        print(e.name, e.shape(), e.uri)
    for e in soma.obsp:
        print(e.name, e.shape(), e.uri)
    for e in soma.varp:
        print(e.name, e.shape(), e.uri)


def test_not_exists():
    soma = tiledbsc.SOMA("/nonesuch/nowhere/never", verbose=False)
    assert not soma.exists()
    assert not soma.obs.exists()
