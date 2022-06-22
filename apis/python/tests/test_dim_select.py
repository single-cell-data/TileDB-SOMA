import tiledbsc
import tiledbsc.io
import tiledb
import anndata
import pandas as pd
import numpy as np

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
    tiledbsc.io.from_anndata(soma, adata)

    assert soma.obs.ids() == [
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

    assert soma.var.ids() == [
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
    assert df.shape == (2, 2)

    df = soma.obsm["X_pca"].dim_select([b"AAGCGACTTTGACG", b"AATGCGTGGACGGA"])
    assert df.shape == (2, 19)

    assert soma.X["data"].dim_select([b"AAGCGACTTTGACG"], [b"AKR1C3"]).shape == (1, 1)
    assert soma.X["data"].dim_select(None, [b"AKR1C3"]).shape == (80, 1)
    assert soma.X["data"].dim_select([b"AAGCGACTTTGACG"], None).shape == (20, 1)
    assert soma.X["data"].dim_select(None, None).shape == (1600, 1)

    tempdir.cleanup()


def test_zeroes_handling():

    # Populate data

    obs_ids = [
        "AAATTCGAATCACG",
        "AATGTTGACAGTCA",
        "AGAGATGATCTCGC",
        "CATGGCCTGTGCAT",
        "CCCAACTGCAATCG",
        "CTAACGGAACCGAT",
        "GAACCTGATGAACC",
        "GCGTAAACACGGTT",
        "TTACCATGAATCGC",
        "TTACGTACGTTCAG",
    ]

    var_ids = [
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
        "S100A9",
        "SDPR",
        "TREML1",
    ]

    n_obs = len(obs_ids)
    n_var = len(var_ids)

    cell_types = ["blööd" if obs_id[1] == "A" else "lung" for obs_id in obs_ids]
    feature_names = [
        "ENSG00000999999" if var_id[1] < "M" else "ENSG00000123456"
        for var_id in var_ids
    ]

    # AnnData requires string indices for obs/var
    obs = pd.DataFrame(
        data={
            "obs_id": np.asarray(obs_ids),
            "cell_type": np.asarray(cell_types),
        },
        index=np.arange(n_obs).astype(str),
    )
    obs.set_index("obs_id", inplace=True)
    var = pd.DataFrame(
        data={
            "var_id": np.asarray(var_ids),
            "feature_name": np.asarray(feature_names),
        },
        index=np.arange(n_var).astype(str),
    )
    var.set_index("var_id", inplace=True)

    X = np.zeros((n_obs, n_var))
    X[0, 0] = 1

    ann = anndata.AnnData(X=X, obs=obs, var=var, dtype=X.dtype)

    # Write SOMA
    tempdir = tempfile.TemporaryDirectory()
    soma_path = tempdir.name
    soma = tiledbsc.SOMA(soma_path)
    tiledbsc.io.from_anndata(soma, ann)

    assert soma.obs.df().shape == (10, 1)
    assert soma.var.df().shape == (16, 1)
    assert soma.X.data.df().shape == (1, 1)  # sparse representation -- IJV triples
    assert soma.X.data.csr().shape == (10, 16)

    tempdir.cleanup()
