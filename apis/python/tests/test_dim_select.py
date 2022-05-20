import anndata
import tiledb
import tiledbsc

import pytest
import tempfile
import os
from pathlib import Path

HERE = Path(__file__).parent


@pytest.fixture
def h5ad_file(request):
    # Tests in this file rely on specific values form this particular input data file.
    input_path = HERE.parent / "anndata/pbmc-small.h5ad"
    return input_path


@pytest.fixture
def adata(h5ad_file):
    return anndata.read_h5ad(h5ad_file)


def test_dim_select(adata):

    # Set up anndata input path and tiledb-group output path
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name

    # Ingest
    soma = tiledbsc.SOMA(output_path, verbose=True)
    soma.from_anndata(adata)

    assert soma.obs.ids() == [
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

    assert soma.var.ids() == [
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

    df = soma.obs.dim_select([b"AAGCGACTTTGACG", b"AATGCGTGGACGGA"])
    assert df.shape == (2, 7)
    assert df.at["AAGCGACTTTGACG", "groups"] == "g1"
    assert df.at["AATGCGTGGACGGA", "nFeature_RNA"] == 73
    #                 orig.ident  nCount_RNA  nFeature_RNA  RNA_snn_res.0.8  letter.idents groups  RNA_snn_res.1
    # obs_id
    # AAGCGACTTTGACG           0       443.0            77                1              1     g1              1
    # AATGCGTGGACGGA           0       389.0            73                1              1     g1              1
    assert soma.obs.dim_select(None).shape == (80, 7)

    df = soma.var.dim_select([b"AKR1C3", b"MYL9"])
    assert df.shape == (2, 5)
    assert df.at["AKR1C3", "vst.variable"] == 1
    assert df.at["MYL9", "vst.variable"] == 1
    assert soma.var.dim_select(None).shape == (20, 5)

    assert sorted(soma.obsm.keys()) == sorted(["X_tsne", "X_pca"])

    df = soma.obsm["X_tsne"].dim_select([b"AAGCGACTTTGACG", b"AATGCGTGGACGGA"])
    assert df.shape == (2, 3)

    df = soma.obsm["X_pca"].dim_select([b"AAGCGACTTTGACG", b"AATGCGTGGACGGA"])
    assert df.shape == (2, 20)

    assert soma.X.data.dim_select([b"AAGCGACTTTGACG"], [b"AKR1C3"]).shape == (1, 3)
    assert soma.X.data.dim_select(None, [b"AKR1C3"]).shape == (80, 3)
    assert soma.X.data.dim_select([b"AAGCGACTTTGACG"], None).shape == (20, 3)
    assert soma.X.data.dim_select(None, None).shape == (1600, 3)

    tempdir.cleanup()
