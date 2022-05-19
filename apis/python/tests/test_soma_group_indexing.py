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

    assert set(soma.get_member_names()) == set(
        ["uns", "varm", "X", "raw", "obsp", "varp", "var", "obsm", "obs"]
    )
    assert set(soma.X.get_member_names()) == set(["data"])
    assert soma.X.data.get_dim_names() == ["obs_id", "var_id"]

    assert soma.obs.get_dim_names() == ["obs_id"]
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

    assert soma.var.get_dim_names() == ["var_id"]
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

    assert set(soma.obsm.get_member_names()) == set(["X_pca", "X_tsne"])
    assert isinstance(soma.obsm["X_pca"], tiledbsc.AnnotationMatrix)
    assert soma.obsm["nonesuch"] is None
    assert soma.obsm["X_pca"].get_dim_names() == ["obs_id"]

    assert set(soma.varm.get_member_names()) == set(["PCs"])
    assert isinstance(soma.varm["PCs"], tiledbsc.AnnotationMatrix)
    assert soma.varm["nonesuch"] is None
    assert soma.varm.get_member_names() == ["PCs"]
    assert soma.varm["PCs"].get_dim_names() == ["var_id"]

    assert set(soma.obsp.get_member_names()) == set(["distances"])
    assert isinstance(soma.obsp["distances"], tiledbsc.AnnotationPairwiseMatrix)
    assert soma.varp["nonesuch"] is None
    assert soma.obsp["distances"].get_dim_names() == ["obs_id_i", "obs_id_j"]

    assert set(soma.uns.get_member_names()) == set(["neighbors"])
    assert isinstance(soma.uns["neighbors"], tiledbsc.UnsGroup)
    assert set(soma.uns["neighbors"].get_member_names()) == set(["params"])
    assert isinstance(soma.uns["neighbors"]["params"], tiledbsc.UnsGroup)
    assert set(soma.uns["neighbors"]["params"].get_member_names()) == set(["method"])
    assert isinstance(soma.uns["neighbors"]["params"]["method"], tiledbsc.UnsArray)
    assert soma.uns["nonesuch"] is None
