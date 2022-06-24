import tempfile
from pathlib import Path

import anndata
import pytest

import tiledbsc
import tiledbsc.io

HERE = Path(__file__).parent


@pytest.fixture
def h5ad_file(request):
    # Tests in this file rely on specific values form this particular input data file.
    input_path = HERE.parent / "anndata/pbmc-small.h5ad"
    return input_path


@pytest.fixture
def adata(h5ad_file):
    return anndata.read_h5ad(h5ad_file)


def test_query(adata):

    # Set up anndata input path and tiledb-group output path
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name

    # Ingest
    soma = tiledbsc.SOMA(output_path)
    tiledbsc.io.from_anndata(soma, adata)

    output = soma.obs.query(
        query_string="nCount_RNA > 10",
        attrs=["nCount_RNA", "orig.ident", "nFeature_RNA"],
    )
    assert output.shape == (80, 3)
    assert output.at["TTACGTACGTTCAG", "nFeature_RNA"] == 39
    #                 orig.ident  nFeature_RNA
    # obs_id
    # AAATTCGAATCACG           0            62
    # AAGCAAGAGCTTAG           0            48
    # AAGCGACTTTGACG           0            77
    # AATGCGTGGACGGA           0            73
    # AATGTTGACAGTCA           0            41
    # ...                    ...           ...
    # TTACGTACGTTCAG           0            39
    # TTGAGGACTACGCA           0            88
    # TTGCATTGAGCTAC           0            40
    # TTGGTACTGAATCC           0            45
    # TTTAGCTGTACTCT           0            86

    # Note: attribute names with "." in them need to be written like 'attr("this.has.dots")' within
    # query-condition expressions.
    output = soma.var.query(
        query_string='attr("vst.mean") > 1',
        attrs=["vst.mean", "vst.variance", "vst.variable"],
    )
    assert output.shape == (9, 3)
    assert output.at["S100A8", "vst.variable"] == 1
    #           vst.mean  vst.variance  vst.variance.expected  vst.variance.standardized  vst.variable
    # var_id
    # GNLY        2.4000     43.762025              24.078566                   1.817468             1
    # HLA-DPB1    7.6500    309.800000             199.154430                   1.555577             1
    # HLA-DQA1    1.9875     32.164399              19.810998                   1.623563             1
    # PF4         2.1500     65.243038              21.482754                   1.939028             1
    # PPBP        5.3250    231.817089              93.148559                   2.488681             1
    # S100A8      2.9750     53.012025              32.261112                   1.643218             1
    # S100A9      6.8375    240.644146             156.219313                   1.540425             1
    # SDPR        1.1125     15.746677               8.835280                   1.680686             1
    # VDAC3       1.1250     30.971519               8.986513                   2.137607             1

    tempdir.cleanup()
