from __future__ import annotations

import json
import random
import tempfile
from copy import deepcopy
from pathlib import Path
from typing import Tuple

import anndata
import numpy as np
import pandas as pd
import pyarrow as pa
import pytest
import scipy
import somacore
from anndata import AnnData
from scipy.sparse import csr_matrix

import tiledbsoma
import tiledbsoma.io
from tiledbsoma import Experiment, _constants, _factory
from tiledbsoma._soma_object import SOMAObject
from tiledbsoma.io._common import _TILEDBSOMA_TYPE, UnsDict, UnsMapping

from ._util import TESTDATA, assert_adata_equal, make_pd_df


@pytest.fixture
def h5ad_file_with_obsm_holes(request):
    # This has zeroes in an obsm matrix so nnz is not num_rows * num_cols
    return TESTDATA / "pbmc3k-with-obsm-zero.h5ad"


@pytest.fixture
def h5ad_file_uns_string_arrays(request):
    # This has uns["louvain_colors"] with dtype.char == "U".
    # It also has uns["more_colors"] in the form '[[...]]', as often occurs in the wild.
    return TESTDATA / "pbmc3k.h5ad"


@pytest.fixture
def h5ad_file_categorical_int_nan(request):
    # This has obs["categ_int_nan"] as a categorical int but with math.nan as a
    # "not-in-the-category" indicator. Such H5AD files do arise in the wild.
    #
    # Reference:
    #   import anndata as ad
    #   import pandas  as pd
    #   import math
    #   adata = adata.read_h5ad("whatever.h5ad")
    #   s = pd.Series(list(range(80)), dtype="category")
    #   s[0] = math.nan
    #   adata.obs["categ_int_nan"] = s
    #   adata.write_h5ad("categorical_int_nan.h5ad")
    return TESTDATA / "categorical_int_nan.h5ad"


@pytest.fixture
def h5ad_file_X_empty(request):
    """adata.X is a zero-cell sparse matrix"""
    return TESTDATA / "x-empty.h5ad"


@pytest.fixture
def h5ad_file_X_none(request):
    """
    adata.X has Python value None if read in non-backed mode; if read in backed
    mode, adata.X is not present as an attribute of adata.
    """
    return TESTDATA / "x-none.h5ad"


@pytest.mark.parametrize(
    "ingest_modes",
    [
        # Standard ingest: normal use-case:
        ["write"],
        # Schema only:
        ["schema_only"],
        # Schema only, then populate:
        ["schema_only", "resume"],
        # User writes data, then a subsequent write creates nothing new:
        ["write", "resume"],
        # "Resume" after no write at all does write new data:
        ["resume"],
    ],
)
@pytest.mark.parametrize(
    "X_kind",
    [tiledbsoma.SparseNDArray, tiledbsoma.DenseNDArray],
)
def test_import_anndata(conftest_pbmc_small, ingest_modes, X_kind):
    original = conftest_pbmc_small.copy()
    conftest_pbmc_small = conftest_pbmc_small.copy()

    have_ingested = False

    tempdir = tempfile.TemporaryDirectory(prefix="test_import_anndata_")
    output_path = tempdir.name

    conftest_pbmc_small.layers["plus1"] = conftest_pbmc_small.X + 1
    orig = conftest_pbmc_small.copy()

    metakey = _constants.SOMA_OBJECT_TYPE_METADATA_KEY  # keystroke-saver
    all2d = (slice(None), slice(None))  # keystroke-saver

    for ingest_mode in ingest_modes:
        uri = tiledbsoma.io.from_anndata(
            output_path,
            orig,
            "RNA",
            ingest_mode=ingest_mode,
            X_kind=X_kind,
        )
        if ingest_mode != "schema_only":
            have_ingested = True

        assert_adata_equal(original, conftest_pbmc_small)

        exp = tiledbsoma.Experiment.open(uri)

        assert exp.metadata[metakey] == "SOMAExperiment"

        # Check obs
        obs = exp.obs.read().concat().to_pandas()
        assert sorted(obs.columns.to_list()) == sorted(
            orig.obs_keys() + ["soma_joinid", "obs_id"]
        )
        assert exp.obs.metadata.get(metakey) == "SOMADataFrame"
        if have_ingested:
            assert sorted(obs["obs_id"]) == sorted(list(orig.obs_names))
        else:
            assert sorted(obs["obs_id"]) == []
        # Convenience accessor
        assert sorted(exp.obs.keys()) == sorted(
            list(orig.obs.keys()) + ["soma_joinid", "obs_id"]
        )

        # Check ms
        assert exp.ms.metadata.get(metakey) == "SOMACollection"
        assert exp.ms["RNA"].metadata.get(metakey) == "SOMAMeasurement"

        # Check var
        var = exp.ms["RNA"].var.read().concat().to_pandas()
        assert sorted(var.columns.to_list()) == sorted(
            orig.var_keys() + ["soma_joinid", "var_id"]
        )
        assert exp.ms["RNA"].var.metadata.get(metakey) == "SOMADataFrame"
        if have_ingested:
            assert sorted(var["var_id"]) == sorted(list(orig.var_names))
        else:
            assert sorted(var["var_id"]) == []
        # Convenience accessor
        assert sorted(exp.ms["RNA"].var.keys()) == sorted(
            list(orig.var.keys()) + ["soma_joinid", "var_id"]
        )

        # Check Xs
        assert exp.ms["RNA"].X.metadata.get(metakey) == "SOMACollection"

        # Check Xs
        assert exp.ms["RNA"].X.metadata.get(metakey) == "SOMACollection"

        # Check X/data (dense in the H5AD)
        for key in ["data", "plus1"]:
            X = exp.ms["RNA"].X[key]
            if X_kind == tiledbsoma.SparseNDArray:
                assert X.metadata.get(metakey) == "SOMASparseNDArray"
                table = X.read(coords=all2d).tables().concat()
                if have_ingested:
                    assert table.shape[0] == orig.X.shape[0] * orig.X.shape[1]
                else:
                    assert table.shape[0] == 0
            else:
                assert X.metadata.get(metakey) == "SOMADenseNDArray"
                matrix = X.read(coords=all2d)
                assert matrix.size == orig.X.size

        # Check raw/X/data (sparse)
        assert exp.ms["raw"].X["data"].metadata.get(metakey) == "SOMASparseNDArray"
        if have_ingested:
            table = exp.ms["raw"].X["data"].read(coords=all2d).tables().concat()
            assert table.shape[0] == orig.raw.X.nnz

        # Check some annotation matrices

        obsm = exp.ms["RNA"].obsm
        assert sorted(obsm.keys()) == sorted(orig.obsm.keys())
        for key, orig_value in orig.obsm.items():
            assert isinstance(obsm[key], somacore.NDArray)
            if have_ingested:
                matrix = obsm[key].read(coords=all2d)
                if not obsm[key].is_sparse:
                    assert matrix.shape == orig_value.shape
            else:
                if isinstance(obsm[key], tiledbsoma.DenseNDArray):
                    with pytest.raises(ValueError):
                        matrix = obsm[key].read(coords=all2d)

        varm = exp.ms["RNA"].varm
        assert sorted(varm.keys()) == sorted(orig.varm.keys())
        for key, orig_value in orig.varm.items():
            assert isinstance(varm[key], somacore.NDArray)
            if have_ingested:
                matrix = varm[key].read(coords=all2d)
                if not varm[key].is_sparse:
                    assert matrix.shape == orig_value.shape
            else:
                if not varm[key].is_sparse:
                    with pytest.raises(ValueError):
                        matrix = varm[key].read(coords=all2d)

        obsp = exp.ms["RNA"].obsp
        assert sorted(obsp.keys()) == sorted(orig.obsp.keys())
        for key, orig_value in orig.obsp.items():
            assert isinstance(obsp[key], somacore.NDArray)
            if have_ingested:
                table = obsp[key].read(coords=all2d).tables().concat()
                assert table.shape[0] == orig_value.nnz

        # pbmc-small has no varp

        tempdir.cleanup()


@pytest.mark.parametrize(
    "X_layer_name",
    [
        None,
        "data",
        "othername",
    ],
)
def test_named_X_layers(conftest_pbmc_small_h5ad_path, X_layer_name):
    tempdir = tempfile.TemporaryDirectory(prefix="test_named_X_layers_")
    soma_path = tempdir.name

    if X_layer_name is None:
        tiledbsoma.io.from_h5ad(
            soma_path,
            conftest_pbmc_small_h5ad_path.as_posix(),
            "RNA",
            ingest_mode="write",
        )
    else:
        tiledbsoma.io.from_h5ad(
            soma_path,
            conftest_pbmc_small_h5ad_path.as_posix(),
            "RNA",
            ingest_mode="write",
            X_layer_name=X_layer_name,
            raw_X_layer_name=X_layer_name,
        )

    with tiledbsoma.Experiment.open(soma_path) as exp:
        if X_layer_name is None:
            assert "data" in exp.ms["RNA"].X
            assert "data" in exp.ms["raw"].X
        else:
            assert X_layer_name in exp.ms["RNA"].X
            assert X_layer_name in exp.ms["raw"].X


def _get_fragment_count(array_uri):
    fragment_uri = Path(array_uri) / "__fragments"
    return len(list(fragment_uri.iterdir())) if fragment_uri.exists() else 0


@pytest.mark.parametrize(
    "resume_mode_h5ad_file",
    [
        TESTDATA / "pbmc-small-x-dense.h5ad",
        TESTDATA / "pbmc-small-x-csr.h5ad",
        TESTDATA / "pbmc-small-x-csc.h5ad",
    ],
)
def test_resume_mode(resume_mode_h5ad_file):
    """
    Makes sure resume-mode ingest after successful ingest of the same input data does not write
    anything new
    """

    tempdir1 = tempfile.TemporaryDirectory(prefix="test_resume_mode_1_")
    output_path1 = tempdir1.name
    tiledbsoma.io.from_h5ad(
        output_path1, resume_mode_h5ad_file.as_posix(), "RNA", ingest_mode="write"
    )

    tempdir2 = tempfile.TemporaryDirectory(prefix="test_resume_mode_2_")
    output_path2 = tempdir2.name
    tiledbsoma.io.from_h5ad(
        output_path2, resume_mode_h5ad_file.as_posix(), "RNA", ingest_mode="write"
    )
    tiledbsoma.io.from_h5ad(
        output_path2, resume_mode_h5ad_file.as_posix(), "RNA", ingest_mode="resume"
    )

    exp1 = _factory.open(output_path1)
    exp2 = _factory.open(output_path2)
    with exp1, exp2:
        assert _get_fragment_count(exp1.obs.uri) == _get_fragment_count(exp2.obs.uri)
        assert _get_fragment_count(exp1.ms["RNA"].var.uri) == _get_fragment_count(
            exp2.ms["RNA"].var.uri
        )
        assert _get_fragment_count(exp1.ms["RNA"].X["data"].uri) == _get_fragment_count(
            exp2.ms["RNA"].X["data"].uri
        )

        meas1 = exp1.ms["RNA"]
        meas2 = exp2.ms["RNA"]

        if "obsm" in meas1:
            for key in meas1.obsm.keys():
                assert _get_fragment_count(meas1.obsm[key].uri) == _get_fragment_count(
                    meas2.obsm[key].uri
                )
        if "varm" in meas1:
            for key in meas1.varm.keys():
                assert _get_fragment_count(meas1.obsm[key].uri) == _get_fragment_count(
                    meas2.obsm[key].uri
                )
        if "obsp" in meas1:
            for key in meas1.obsp.keys():
                assert _get_fragment_count(meas1.obsp[key].uri) == _get_fragment_count(
                    meas2.obsp[key].uri
                )
        if "varp" in meas1:
            for key in meas1.varp.keys():
                assert _get_fragment_count(meas1.varm[key].uri) == _get_fragment_count(
                    meas2.varm[key].uri
                )


@pytest.mark.parametrize("use_relative_uri", [False, True, None])
def test_ingest_relative(conftest_pbmc3k_h5ad_path, use_relative_uri):
    tempdir = tempfile.TemporaryDirectory(prefix="test_ingest_relative_")
    output_path = tempdir.name

    tiledbsoma.io.from_h5ad(
        output_path,
        conftest_pbmc3k_h5ad_path.as_posix(),
        measurement_name="RNA",
        use_relative_uri=use_relative_uri,
    )

    # * If they ask for relative=True, they should get that.
    # * If they ask for relative=False, they should get that.
    # * If they ask for relative=None, they should get the default which, for local disk (these
    #   tests) is True.
    expected_relative = use_relative_uri
    if use_relative_uri is None:
        expected_relative = True  # since local disk

    with tiledbsoma.Experiment.open(output_path) as G:
        assert G._handle._handle.is_relative("obs") == expected_relative
        assert G._handle._handle.is_relative("ms") == expected_relative

        assert G.ms._handle._handle.is_relative("RNA") == expected_relative
        assert G.ms["RNA"]._handle._handle.is_relative("var") == expected_relative
        assert G.ms["RNA"]._handle._handle.is_relative("X") == expected_relative
        assert G.ms["RNA"].X._handle._handle.is_relative("data") == expected_relative

        for collection_name in [
            "obsm",
            "obsp",
            "varm",
        ]:  # conftest_h5ad_file_extended has no varp
            for member in G.ms["RNA"][collection_name]:
                assert (
                    G.ms["RNA"][collection_name]._handle._handle.is_relative(member)
                    == expected_relative
                )

        assert G.ms._handle._handle.is_relative("raw") == expected_relative
        assert G.ms["raw"]._handle._handle.is_relative("var") == expected_relative
        assert G.ms["raw"]._handle._handle.is_relative("X") == expected_relative
        assert G.ms["raw"].X._handle._handle.is_relative("data") == expected_relative


@pytest.mark.parametrize("ingest_uns_keys", [["louvain_colors"], None])
def test_ingest_uns(
    tmp_path: Path,
    conftest_pbmc3k_h5ad_path,
    conftest_pbmc3k_adata,
    ingest_uns_keys,
):
    tmp_uri = tmp_path.as_uri()
    adata_extended2 = anndata.read_h5ad(conftest_pbmc3k_h5ad_path)
    uri = tiledbsoma.io.from_anndata(
        tmp_uri,
        adata_extended2,
        measurement_name="hello",
        uns_keys=ingest_uns_keys,
    )

    assert_adata_equal(conftest_pbmc3k_adata, adata_extended2)

    with tiledbsoma.Experiment.open(uri) as exp:
        uns = exp.ms["hello"]["uns"]
        assert isinstance(uns, tiledbsoma.Collection)
        assert uns.metadata[_TILEDBSOMA_TYPE] == "uns"
        if ingest_uns_keys is None:
            assert set(uns) == {
                "draw_graph",
                "louvain",
                "louvain_colors",
                "neighbors",
                "pca",
                "rank_genes_groups",
            }
            assert isinstance(uns["louvain_colors"], tiledbsoma.DataFrame)
            rgg = uns["rank_genes_groups"]
            assert set(rgg) == {"params"}, "structured arrays not imported"
            assert rgg["params"].metadata.items() >= {
                ("groupby", "louvain"),
                ("method", "t-test_overestim_var"),
                ("reference", "rest"),
            }
            dg_params = uns["draw_graph"]["params"]
            assert isinstance(dg_params, tiledbsoma.Collection)
            assert dg_params.metadata["layout"] == "fr"
            random_state = dg_params["random_state"]
            assert isinstance(random_state, tiledbsoma.DenseNDArray)
            assert np.array_equal(random_state.read().to_numpy(), np.array([0]))
            got_pca_variance = uns["pca"]["variance"].read().to_numpy()
            assert np.array_equal(
                got_pca_variance, adata_extended2.uns["pca"]["variance"]
            )
        else:
            assert set(uns) == set(ingest_uns_keys)


def test_ingest_uns_string_arrays(h5ad_file_uns_string_arrays):
    tempdir = tempfile.TemporaryDirectory(prefix="test_ingest_uns_string_arrays_")
    output_path = tempdir.name

    tiledbsoma.io.from_h5ad(
        output_path,
        h5ad_file_uns_string_arrays.as_posix(),
        measurement_name="RNA",
    )

    with tiledbsoma.Experiment.open(output_path) as exp:
        with tiledbsoma.DataFrame.open(
            exp.ms["RNA"]["uns"]["louvain_colors"].uri
        ) as df:
            contents = df.read().concat()
            assert contents.shape == (8, 2)
            assert len(contents["values"]) == 8
            assert contents["values"][0].as_py() == "#1f77b4"

        with tiledbsoma.DataFrame.open(exp.ms["RNA"]["uns"]["more_colors"].uri) as df:
            contents = df.read().concat()
            assert contents.shape == (8, 2)
            assert len(contents["values_0"]) == 8
            assert contents["values_0"][0].as_py() == "#1f77b4"


def test_add_matrix_to_collection(conftest_pbmc_small):
    tempdir = tempfile.TemporaryDirectory(prefix="test_add_matrix_to_collection_")
    output_path = tempdir.name

    original = conftest_pbmc_small.copy()

    uri = tiledbsoma.io.from_anndata(
        output_path, conftest_pbmc_small, measurement_name="RNA"
    )

    assert_adata_equal(original, conftest_pbmc_small)

    exp = tiledbsoma.Experiment.open(uri)
    with _factory.open(output_path) as exp_r:
        assert list(exp_r.ms["RNA"].X.keys()) == ["data"]
        with pytest.raises(tiledbsoma.SOMAError):
            tiledbsoma.io.add_X_layer(
                exp, "RNA", "data2", conftest_pbmc_small.X
            )  # not open for read
    with _factory.open(output_path, "w") as exp:
        tiledbsoma.io.add_X_layer(exp, "RNA", "data2", conftest_pbmc_small.X)
    with pytest.raises(tiledbsoma.SOMAError):
        tiledbsoma.io.add_X_layer(exp, "RNA", "data3", conftest_pbmc_small.X)  # closed
    with _factory.open(output_path) as exp_r:
        assert sorted(list(exp_r.ms["RNA"].X.keys())) == ["data", "data2"]

    with _factory.open(output_path, "w") as exp:
        with pytest.raises(KeyError):
            tiledbsoma.io.add_X_layer(exp, "nonesuch", "data3", conftest_pbmc_small.X)

    with _factory.open(output_path) as exp_r:
        assert sorted(list(exp_r.ms["RNA"].obsm.keys())) == sorted(
            list(conftest_pbmc_small.obsm.keys())
        )
    with _factory.open(output_path, "w") as exp:
        tiledbsoma.io.add_matrix_to_collection(
            exp, "RNA", "obsm", "X_pcb", conftest_pbmc_small.obsm["X_pca"]
        )
    with _factory.open(output_path) as exp_r:
        assert sorted(list(exp_r.ms["RNA"].obsm.keys())) == sorted(
            list(conftest_pbmc_small.obsm.keys()) + ["X_pcb"]
        )

    with _factory.open(output_path, "w") as exp:
        with pytest.raises(KeyError):
            tiledbsoma.io.add_matrix_to_collection(
                exp, "nonesuch", "obsm", "X_pcc", conftest_pbmc_small.obsm["X_pca"]
            )

        tiledbsoma.io.add_matrix_to_collection(
            exp, "RNA", "newthing", "X_pcd", conftest_pbmc_small.obsm["X_pca"]
        )
    with _factory.open(output_path) as exp_r:
        assert sorted(list(exp_r.ms["RNA"]["newthing"].keys())) == ["X_pcd"]


# This tests tiledbsoma.io._create_or_open_coll -> tiledbsoma.io._create_or_open_collection
# compatibility. As a rule we do not offer any backward-compatibility assurances for non-public
# modules or methods -- those whose names start with an underscore.  For this single case we are
# making an exception. For future code-imitation purposes, please be aware this is a pattern to be
# avoided in the future, not imitated.
def test_add_matrix_to_collection_1_2_7(conftest_pbmc_small):
    def add_X_layer(
        exp: tiledbsoma.Experiment,
        measurement_name: str,
        X_layer_name: str,
        X_layer_data,  # E.g. a scipy.csr_matrix from scanpy analysis
        ingest_mode: str = "write",
        use_relative_uri: bool | None = None,
    ) -> None:
        if exp.closed or exp.mode != "w":
            raise tiledbsoma.SOMAError(f"Experiment must be open for write: {exp.uri}")
        add_matrix_to_collection(
            exp,
            measurement_name,
            "X",
            X_layer_name,
            X_layer_data,
            use_relative_uri=use_relative_uri,
        )

    def add_matrix_to_collection(
        exp: tiledbsoma.Experiment,
        measurement_name: str,
        collection_name: str,
        matrix_name: str,
        matrix_data,  # E.g. a scipy.csr_matrix from scanpy analysis
        ingest_mode: str = "write",
        use_relative_uri: bool | None = None,
        context: tiledbsoma.SOMATileDBContext | None = None,
    ) -> None:
        # For local disk and S3, creation and storage URIs are identical.  For
        # cloud, creation URIs look like tiledb://namespace/s3://bucket/path/to/obj
        # whereas storage URIs (for the same object) look like
        # tiledb://namespace/uuid.  When the caller passes a creation URI (which
        # they must) via exp.uri, we need to follow that.
        extend_creation_uri = exp.uri.startswith("tiledb://")

        with exp.ms[measurement_name] as meas:
            if extend_creation_uri:
                coll_uri = f"{exp.uri}/ms/{measurement_name}/{collection_name}"
            else:
                coll_uri = f"{meas.uri}/{collection_name}"

            if collection_name in meas:
                coll = meas[collection_name]
            else:
                coll = tiledbsoma.io.ingest._create_or_open_coll(
                    tiledbsoma.Collection,
                    coll_uri,
                    ingest_mode=ingest_mode,
                    context=context,
                )
                tiledbsoma.io.ingest._maybe_set(
                    meas, collection_name, coll, use_relative_uri=use_relative_uri
                )
            with coll:
                matrix_uri = f"{coll_uri}/{matrix_name}"

                with tiledbsoma.io.ingest.create_from_matrix(
                    tiledbsoma.SparseNDArray,
                    matrix_uri,
                    matrix_data,
                    context=context,
                ) as sparse_nd_array:
                    tiledbsoma.io.ingest._maybe_set(
                        coll,
                        matrix_name,
                        sparse_nd_array,
                        use_relative_uri=use_relative_uri,
                    )

    tempdir = tempfile.TemporaryDirectory(prefix="test_add_matrix_to_collection_1_2_7_")
    output_path = tempdir.name
    original = conftest_pbmc_small.copy()

    uri = tiledbsoma.io.from_anndata(
        output_path, conftest_pbmc_small, measurement_name="RNA"
    )

    assert_adata_equal(original, conftest_pbmc_small)

    exp = tiledbsoma.Experiment.open(uri)
    with _factory.open(output_path) as exp_r:
        assert list(exp_r.ms["RNA"].X.keys()) == ["data"]
        with pytest.raises(tiledbsoma.SOMAError):
            add_X_layer(exp, "RNA", "data2", conftest_pbmc_small.X)  # not open for read
    with _factory.open(output_path, "w") as exp:
        add_X_layer(exp, "RNA", "data2", conftest_pbmc_small.X)
    with pytest.raises(tiledbsoma.SOMAError):
        add_X_layer(exp, "RNA", "data3", conftest_pbmc_small.X)  # closed
    with _factory.open(output_path) as exp_r:
        assert sorted(list(exp_r.ms["RNA"].X.keys())) == ["data", "data2"]

    with _factory.open(output_path, "w") as exp:
        with pytest.raises(KeyError):
            add_X_layer(exp, "nonesuch", "data3", conftest_pbmc_small.X)

    with _factory.open(output_path) as exp_r:
        assert sorted(list(exp_r.ms["RNA"].obsm.keys())) == sorted(
            list(conftest_pbmc_small.obsm.keys())
        )

    with _factory.open(output_path, "w") as exp:
        add_matrix_to_collection(
            exp, "RNA", "obsm", "X_pcb", conftest_pbmc_small.obsm["X_pca"]
        )
    with _factory.open(output_path) as exp_r:
        assert sorted(list(exp_r.ms["RNA"].obsm.keys())) == sorted(
            list(conftest_pbmc_small.obsm.keys()) + ["X_pcb"]
        )

    with _factory.open(output_path, "w") as exp:
        # It's nonsense biologically to add this to varp, but as a fake-test unit-test case, we can
        # use varp to test adding to a not-yet-existing collection.
        add_matrix_to_collection(
            exp, "RNA", "varp", "X_pcb", conftest_pbmc_small.obsm["X_pca"]
        )
    with _factory.open(output_path) as exp_r:
        assert sorted(list(exp_r.ms["RNA"].varp.keys())) == sorted(
            list(conftest_pbmc_small.varp.keys()) + ["X_pcb"]
        )

    with _factory.open(output_path, "w") as exp:
        with pytest.raises(KeyError):
            add_matrix_to_collection(
                exp, "nonesuch", "obsm", "X_pcc", conftest_pbmc_small.obsm["X_pca"]
            )

        add_matrix_to_collection(
            exp, "RNA", "newthing", "X_pcd", conftest_pbmc_small.obsm["X_pca"]
        )
    with _factory.open(output_path) as exp_r:
        assert sorted(list(exp_r.ms["RNA"]["newthing"].keys())) == ["X_pcd"]


def test_export_anndata(conftest_pbmc_small):
    tempdir = tempfile.TemporaryDirectory(prefix="test_export_anndata_")
    output_path = tempdir.name

    original = conftest_pbmc_small.copy()

    tiledbsoma.io.from_anndata(output_path, conftest_pbmc_small, measurement_name="RNA")

    assert_adata_equal(original, conftest_pbmc_small)

    with _factory.open(output_path) as exp:
        with pytest.raises(ValueError):
            tiledbsoma.io.to_anndata(
                exp, measurement_name="RNA", obs_id_name="nonesuch"
            )
        with pytest.raises(ValueError):
            tiledbsoma.io.to_anndata(
                exp, measurement_name="RNA", var_id_name="nonesuch"
            )
        with pytest.raises(ValueError):
            tiledbsoma.io.to_anndata(exp, measurement_name="nonesuch")
        with pytest.raises(ValueError):
            tiledbsoma.io.to_anndata(
                exp, measurement_name="RNA", X_layer_name="nonesuch"
            )

        readback = tiledbsoma.io.to_anndata(exp, measurement_name="RNA")

    assert readback.obs.shape == conftest_pbmc_small.obs.shape
    assert readback.var.shape == conftest_pbmc_small.var.shape
    assert readback.X.shape == conftest_pbmc_small.X.shape

    for key in conftest_pbmc_small.obsm.keys():
        assert readback.obsm[key].shape == conftest_pbmc_small.obsm[key].shape
    for key in conftest_pbmc_small.varm.keys():
        assert readback.varm[key].shape == conftest_pbmc_small.varm[key].shape
    for key in conftest_pbmc_small.obsp.keys():
        assert readback.obsp[key].shape == conftest_pbmc_small.obsp[key].shape
    for key in conftest_pbmc_small.varp.keys():
        assert readback.varp[key].shape == conftest_pbmc_small.varp[key].shape

    with _factory.open(output_path) as exp:
        readback = tiledbsoma.io.to_anndata(
            exp, measurement_name="RNA", X_layer_name=None
        )
        assert readback.obs.shape == conftest_pbmc_small.obs.shape
        assert readback.var.shape == conftest_pbmc_small.var.shape
        assert readback.X is None


def test_ingest_additional_metadata(conftest_pbmc_small):
    tempdir = tempfile.TemporaryDirectory(prefix="test_ingest_additional_metadata_")
    output_path = tempdir.name

    additional_metadata = {"key1": "val1", "key2": "val2"}

    tiledbsoma.io.from_anndata(
        output_path,
        conftest_pbmc_small,
        measurement_name="RNA",
        additional_metadata=additional_metadata,
    )

    def check(tdbo: SOMAObject):
        for k, v in additional_metadata.items():
            assert tdbo.metadata[k] == v

    with _factory.open(output_path, soma_type=Experiment) as exp:
        check(exp)
        check(exp.obs)
        rna = exp.ms["RNA"]
        check(rna)
        check(rna.X)
        check(rna.obsm)
        check(rna.obsp)
        check(rna.varm)
        # No varp in pbmc-small.h5ad

        raw = exp.ms["raw"]
        check(raw)
        check(raw.X)


def test_null_obs(conftest_pbmc_small, tmp_path: Path):
    output_path = tmp_path.as_uri()
    seed = 42
    # Create column of all null values

    # Note: categorical all-nan goes through some type-inference hijinx on the
    # Pandas -> Arrow path as detailed in
    # https://github.com/single-cell-data/TileDB-SOMA/issues/3685.
    # The TL;DR is that:
    # * If one is not careful, it type-infers as categorical-of-double.
    # * If one is very careful, it can be made to type-infer as categorical-of-string.
    # * The latter case is by far the most common case, and we must bias toward it.
    conftest_pbmc_small.obs["nan_categorical_all"] = pd.Categorical(
        [np.nan] * conftest_pbmc_small.n_obs,
        dtype=pd.CategoricalDtype(categories=[], ordered=False),
    )
    conftest_pbmc_small.obs["none_categorical_all"] = pd.Categorical(
        [None] * conftest_pbmc_small.n_obs,
        dtype=pd.CategoricalDtype(categories=[], ordered=False),
    )
    conftest_pbmc_small.obs["empty_extension_all"] = pd.Series(
        [np.nan] * conftest_pbmc_small.n_obs, dtype=pd.Int64Dtype()
    )
    # Create column of partially-null values
    rng = np.random.RandomState(seed)

    conftest_pbmc_small.obs["empty_categorical_partial"] = rng.choice(
        (np.nan, 1.0), conftest_pbmc_small.n_obs, True
    )
    conftest_pbmc_small.obs["empty_extension_partial"] = pd.Series(
        [1] * conftest_pbmc_small.n_obs + [np.nan], dtype=pd.Int64Dtype()
    )

    original = conftest_pbmc_small.copy()
    uri = tiledbsoma.io.from_anndata(
        output_path,
        conftest_pbmc_small,
        "RNA",
        ingest_mode="write",
        X_kind=tiledbsoma.SparseNDArray,
    )
    assert_adata_equal(original, conftest_pbmc_small)

    with tiledbsoma.Experiment.open(uri) as exp:
        getter = exp.obs.schema.field

        # Explicitly check columns created above
        assert getter("nan_categorical_all").nullable
        assert getter("none_categorical_all").nullable
        assert getter("empty_categorical_partial").nullable
        assert getter("empty_extension_all").nullable
        assert getter("empty_extension_partial").nullable

        # https://github.com/single-cell-data/TileDB-SOMA/issues/3685.
        assert getter("nan_categorical_all").type.value_type == pa.string()
        assert getter("none_categorical_all").type.value_type == pa.string()

        # For every column in the data frame ensure that `isnullable` reflects
        # he null-ness of the Pandas data frame
        for k in conftest_pbmc_small.obs:
            assert getter(k).nullable


def test_export_obsm_with_holes(h5ad_file_with_obsm_holes, tmp_path):
    adata = anndata.read_h5ad(h5ad_file_with_obsm_holes.as_posix())
    original = adata.copy()

    # This data file is prepared such that obsm["X_pca"] has shape (2638, 50)
    # but its [0][0] element is a 0, so when it's stored as sparse, its nnz
    # is not 2638*50=131900.
    ado = adata.obsm["X_pca"]
    assert ado.shape == (2638, 50)

    output_path = tmp_path.as_posix()
    tiledbsoma.io.from_anndata(output_path, adata, "RNA")

    assert_adata_equal(original, adata)

    with tiledbsoma.Experiment.open(output_path) as exp:
        meta = exp.ms["RNA"].obsm["X_pca"].metadata
        with pytest.raises(KeyError):
            meta["soma_dim_0_domain_lower"]
        with pytest.raises(KeyError):
            meta["soma_dim_0_domain_upper"]
        with pytest.raises(KeyError):
            meta["soma_dim_1_domain_lower"]
        with pytest.raises(KeyError):
            meta["soma_dim_1_domain_upper"]
        assert meta["soma_object_type"] == "SOMASparseNDArray"

        # No longer throws as of new-shape feature in TileDB-SOMA 1.15.
        try3 = tiledbsoma.io.to_anndata(exp, "RNA")
        assert try3.obsm["X_pca"].shape == (2638, 50)

        try4 = tiledbsoma.io.to_anndata(
            exp, "RNA", obsm_varm_width_hints={"obsm": {"X_pca": 50}}
        )
        assert try4.obsm["X_pca"].shape == (2638, 50)


def test_X_empty(h5ad_file_X_empty):
    tempdir = tempfile.TemporaryDirectory(prefix="test_X_empty_")
    output_path = tempdir.name
    tiledbsoma.io.from_h5ad(
        output_path, h5ad_file_X_empty.as_posix(), measurement_name="RNA"
    )

    with tiledbsoma.Experiment.open(output_path) as exp:
        assert exp.obs.count == 2638
        assert exp.ms["RNA"].var.count == 1838
        assert "data" in exp.ms["RNA"].X
        assert exp.ms["RNA"].X["data"].nnz == 0

        tiledbsoma.io.to_anndata(exp, measurement_name="RNA")
        # TODO: more


def test_X_none(h5ad_file_X_none):
    tempdir = tempfile.TemporaryDirectory(prefix="test_X_none_")
    output_path = tempdir.name
    tiledbsoma.io.from_h5ad(
        output_path, h5ad_file_X_none.as_posix(), measurement_name="RNA"
    )

    with tiledbsoma.Experiment.open(output_path) as exp:
        assert exp.obs.count == 2638
        assert exp.ms["RNA"].var.count == 1838
        assert list(exp.ms["RNA"].X.keys()) == []

        adata = tiledbsoma.io.to_anndata(exp, measurement_name="RNA")
        assert adata.obs is not None
        assert adata.var is not None
        assert adata.X is None

        adata = tiledbsoma.io.to_anndata(exp, measurement_name="RNA", X_layer_name=None)
        assert adata.obs is not None
        assert adata.var is not None
        assert adata.X is None

        with pytest.raises(ValueError):
            tiledbsoma.io.to_anndata(exp, measurement_name="RNA", X_layer_name="data")

        with pytest.raises(ValueError):
            adata = tiledbsoma.io.to_anndata(
                exp, measurement_name="RNA", X_layer_name="nonesuch"
            )


# There exist in the wild AnnData files with categorical-int columns where "not in the category" is
# indicated by the presence of floating-point math.NaN in cells. Here we test that we can ingest
# this.
def test_obs_with_categorical_int_nan_enumeration(
    tmp_path, h5ad_file_categorical_int_nan
):
    output_path = tmp_path.as_uri()

    tiledbsoma.io.from_h5ad(
        output_path, h5ad_file_categorical_int_nan, measurement_name="RNA"
    )


@pytest.mark.parametrize("obs_id_name", ["obs_id", "cells_are_great"])
@pytest.mark.parametrize("var_id_name", ["var_id", "genes_are_nice_too"])
@pytest.mark.parametrize("indexify_obs", [True, False])
@pytest.mark.parametrize("indexify_var", [True, False])
def test_id_names(tmp_path, obs_id_name, var_id_name, indexify_obs, indexify_var):
    obs_ids = ["AAAT", "CATG", "CTGA", "TCTG", "TGAG", "TTTG"]
    var_ids = ["AKT1", "APOE", "ESR1", "TP53", "VEGFA", "ZZZ3"]

    n_obs = len(obs_ids)
    n_var = len(var_ids)

    obs = pd.DataFrame(
        data={
            obs_id_name: np.asarray(obs_ids),
            "cell_type": pd.Categorical(
                [["B cell", "T cell"][e % 2] for e in range(n_obs)],
                categories=["B cell", "T cell"],
                ordered=True,
            ),
        },
        index=np.arange(n_obs).astype(str),
    )
    if indexify_obs:
        obs.set_index(obs_id_name, inplace=True)

    var = pd.DataFrame(
        data={
            var_id_name: np.asarray(var_ids),
            "counter": np.asarray(range(n_var), dtype=np.float32),
        },
        index=np.arange(n_var).astype(str),
    )
    if indexify_var:
        var.set_index(var_id_name, inplace=True)

    X = np.zeros([n_obs, n_var])
    for i in range(n_obs):
        for j in range(n_var):
            if (i + j) % 2 == 1:
                X[i, j] = 100 + 10 * i + j

    adata = anndata.AnnData(X=X, obs=obs, var=var, dtype=X.dtype)
    original = adata.copy()

    uri = tmp_path.as_posix()

    # Implicitly, a check for no-throw
    tiledbsoma.io.from_anndata(
        uri,
        adata,
        measurement_name="RNA",
        obs_id_name=obs_id_name,
        var_id_name=var_id_name,
    )
    assert_adata_equal(original, adata)

    with tiledbsoma.Experiment.open(uri) as exp:
        assert obs_id_name in exp.obs.keys()
        assert var_id_name in exp.ms["RNA"].var.keys()

        if indexify_obs:
            expected_obs_keys = ["soma_joinid", obs_id_name] + adata.obs_keys()
        else:
            expected_obs_keys = ["soma_joinid"] + adata.obs_keys()
        assert list(exp.obs.keys()) == expected_obs_keys

        if indexify_var:
            expected_var_keys = ["soma_joinid", var_id_name] + adata.var_keys()
        else:
            expected_var_keys = ["soma_joinid"] + adata.var_keys()
        assert list(exp.ms["RNA"].var.keys()) == expected_var_keys

        # Implicitly, a check for no-throw
        bdata = tiledbsoma.io.to_anndata(
            exp,
            measurement_name="RNA",
            obs_id_name=obs_id_name,
            var_id_name=var_id_name,
        )

        soma_obs = exp.obs.read(column_names=[obs_id_name]).concat().to_pandas()
        soma_var = (
            exp.ms["RNA"].var.read(column_names=[var_id_name]).concat().to_pandas()
        )
        assert list(bdata.obs.index) == list(soma_obs[obs_id_name])
        assert list(bdata.var.index) == list(soma_var[var_id_name])


TEST_UNS: UnsDict = {
    # Scalars are stored in SOMA as metadata
    "int_scalar": 7,
    "float_scalar": 8.5,
    "string_scalar": "hello",
    # pd.DataFrames are stored in SOMA as `DataFrame`s
    "pd_df_indexed": make_pd_df("0,1,2", column_1="d,e,f"),
    "pd_df_nonindexed": make_pd_df(column_1="g,h,i"),
    # `np.ndarray`s are stored in SOMA as `DenseArray`s
    "np_ndarray_1d": np.asarray([1, 2, 3]),
    "np_ndarray_2d": np.asarray([[1, 2, 3], [4, 5, 6]]),
    # Nested dicts are stored in SOMA as `SOMACollection`s
    "strings": {
        # `np.ndarray`s of strings are stored in SOMA as `DataFrame`s, since SOMA ND arrays are
        # necessarily arrays *of numbers*. This is okay since the one and only job of SOMA uns is
        # to faithfully ingest from AnnData and outgest back.
        "string_np_ndarray_1d": np.asarray(["j", "k", "l"]),
        "string_np_ndarray_2d": np.asarray([["m", "n", "o"], ["p", "q", "r"]]),
    },
}


def make_uns_adata(
    tmp_path: Path,
    measurement_name: str = "RNA",
    uns: UnsMapping | None = None,
) -> Tuple[str, AnnData]:
    obs = pd.DataFrame(
        data={"obs_id": np.asarray(["a", "b", "c"])},
        index=np.arange(3).astype(str),
    )
    var = pd.DataFrame(
        data={"var_id": np.asarray(["x", "y"])},
        index=np.arange(2).astype(str),
    )
    X = csr_matrix(np.zeros([3, 2]))

    if uns is None:
        uns = TEST_UNS

    adata = anndata.AnnData(
        obs=obs,
        var=var,
        X=X,
        uns=uns,
        dtype=X.dtype,
    )
    adata0 = deepcopy(adata)

    soma_uri = tmp_path.as_posix()

    tiledbsoma.io.from_anndata(soma_uri, adata, measurement_name=measurement_name)
    assert_adata_equal(adata0, adata)

    return soma_uri, adata


@pytest.mark.parametrize(
    "outgest_uns_keys", [["int_scalar", "strings", "np_ndarray_2d"], None]
)
def test_uns_io(tmp_path, outgest_uns_keys):
    soma_uri, adata = make_uns_adata(tmp_path)

    with tiledbsoma.Experiment.open(soma_uri) as exp:
        adata2 = tiledbsoma.io.to_anndata(
            exp,
            measurement_name="RNA",
            uns_keys=outgest_uns_keys,
        )

    expected_adata = deepcopy(adata)

    # Outgest also fails to restore `obs` and `var` correctly, in this case because the ingested
    # `obs`/`var` had columns named "obs_id"/"var_id", which get mistaken for "default index"
    # columns and set as `df.index` on outgest; their names are also removed. This corresponds to
    # case #2 from https://github.com/single-cell-data/TileDB-SOMA/issues/2829.
    # TODO: fix `to_anndata` to restore `obs` and `var` as ingested.
    expected_adata.obs = expected_adata.obs.set_index("obs_id")
    expected_adata.obs.index.name = None
    expected_adata.var = expected_adata.var.set_index("var_id")
    expected_adata.var.index.name = None

    if outgest_uns_keys is not None:
        expected_adata.uns = {
            k: v for k, v in expected_adata.uns.items() if k in outgest_uns_keys
        }

    assert_adata_equal(expected_adata, adata2)


@pytest.mark.parametrize("write_index", [0, 1])
def test_string_nan_columns(tmp_path, conftest_pbmc_small, write_index):
    # Use case:
    #
    # 1. Anndata has column filled with np.nan
    # 2. Import AnnData to SOMA
    # 3. Export SOMA back to AnnData
    # 4. Fill all/part of empty column with string values

    # Step 1
    conftest_pbmc_small.obs["new_col"] = pd.Series(data=np.nan, dtype=np.dtype(str))

    # Step 2
    uri = tmp_path.as_posix()
    original = conftest_pbmc_small.copy()
    tiledbsoma.io.from_anndata(uri, conftest_pbmc_small, measurement_name="RNA")
    assert_adata_equal(original, conftest_pbmc_small)

    # Step 3
    with tiledbsoma.open(uri, "r") as exp:
        bdata = tiledbsoma.io.to_anndata(exp, measurement_name="RNA")

    # Step 4
    bdata.obs["new_col"][write_index] = "abc"
    with tiledbsoma.open(uri, "w") as exp:
        # Implicit assert here that nothing throws
        tiledbsoma.io.update_obs(exp, bdata.obs)

    # TODO: asserts


@pytest.mark.parametrize("obs_index_name", [None, "obs_id", "cell_id"])
@pytest.mark.parametrize("var_index_name", [None, "var_id", "gene_id"])
def test_index_names_io(tmp_path, obs_index_name, var_index_name):
    nobs = 200
    nvar = 100
    xocc = 0.3
    measurement_name = "meas"

    # White-box-test this, which we leverage inside tiledbsoma.io
    assert json.loads("null") is None

    obs_ids = ["cell_%08d" % (i) for i in range(nobs)]
    var_ids = ["gene_%08d" % (j) for j in range(nvar)]

    cell_types = [["B cell", "T cell"][e % 2] for e in range(nobs)]
    obs = pd.DataFrame(
        data={
            obs_index_name: np.asarray(obs_ids),
            "cell_type": pd.Categorical(cell_types),
        },
        index=np.arange(nobs).astype(str),
    )
    if obs_index_name is not None:
        obs.set_index(obs_index_name, inplace=True)

    var = pd.DataFrame(
        data={
            var_index_name: np.asarray(var_ids),
            "squares": np.asarray([i**2 for i in range(nvar)]),
        },
        index=np.arange(len(var_ids)).astype(str),
    )
    if var_index_name is not None:
        var.set_index(var_index_name, inplace=True)

    X = scipy.sparse.random(nobs, nvar, density=xocc, dtype=np.float64).tocsr()

    adata = anndata.AnnData(X=X, obs=obs, var=var)

    soma_uri = tmp_path.as_posix()

    original = adata.copy()
    tiledbsoma.io.from_anndata(soma_uri, adata, measurement_name)
    assert_adata_equal(original, adata)

    with tiledbsoma.Experiment.open(soma_uri) as exp:
        bdata = tiledbsoma.io.to_anndata(exp, measurement_name)

    if obs_index_name is None:
        assert adata.obs.index.name is None
        assert bdata.obs.index.name is None
    else:
        assert adata.obs.index.name == bdata.obs.index.name

    if var_index_name is None:
        assert adata.var.index.name is None
        assert bdata.var.index.name is None
    else:
        assert adata.var.index.name == bdata.var.index.name


def test_obsm_data_type(conftest_pbmc_small):
    tempdir = tempfile.TemporaryDirectory(prefix="test_obsm_data_type_")
    soma_path = tempdir.name
    bdata = anndata.AnnData(
        X=conftest_pbmc_small.X,
        obs=conftest_pbmc_small.obs,
        var=conftest_pbmc_small.var,
        obsm={"testing": conftest_pbmc_small.obs},
    )

    with pytest.raises(TypeError):
        tiledbsoma.io.from_anndata(soma_path, bdata, measurement_name="RNA")

    assert not any(Path(soma_path).iterdir())


def test_outgest_X_layers(tmp_path):
    nobs = 200
    nvar = 100
    xocc = 0.3
    measurement_name = "meas"

    obs_ids = ["cell_%08d" % (i) for i in range(nobs)]
    var_ids = ["gene_%08d" % (j) for j in range(nvar)]

    cell_types = [["B cell", "T cell"][e % 2] for e in range(nobs)]
    obs = pd.DataFrame(
        data={
            "obs_id": np.asarray(obs_ids),
            "cell_type": pd.Categorical(cell_types),
        },
        index=np.arange(nobs).astype(str),
    )
    obs.set_index("obs_id", inplace=True)

    var = pd.DataFrame(
        data={
            "var_id": np.asarray(var_ids),
            "squares": np.asarray([i**2 for i in range(nvar)]),
        },
        index=np.arange(len(var_ids)).astype(str),
    )
    var.set_index("var_id", inplace=True)

    X = scipy.sparse.random(nobs, nvar, density=xocc, dtype=np.float64).tocsr()
    layers = {"data2": X, "data3": X}

    adata = anndata.AnnData(X=X, obs=obs, var=var, layers=layers)

    soma_uri = tmp_path.as_posix()

    tiledbsoma.io.from_anndata(soma_uri, adata, measurement_name)

    # @pytest.mark.parametrize is intentionally not used here. For this
    # particular test-case, it's simpler to spell out the checks.

    with tiledbsoma.Experiment.open(soma_uri) as exp:
        measurement = exp.ms[measurement_name]

        bdata = tiledbsoma.io.to_anndata(exp, measurement_name)
        assert bdata.X is not None
        assert len(bdata.layers) == 0

        bdata = tiledbsoma.io.to_anndata(
            exp, measurement_name, X_layer_name=None, extra_X_layer_names=[]
        )
        assert bdata.X is None
        assert len(bdata.layers) == 0

        with pytest.raises(ValueError):
            tiledbsoma.io.to_anndata(
                exp, measurement_name, X_layer_name=None, extra_X_layer_names=["data"]
            )

        with pytest.raises(ValueError):
            tiledbsoma.io.to_anndata(exp, measurement_name, X_layer_name="nonesuch")

        with pytest.raises(ValueError):
            tiledbsoma.io.to_anndata(
                exp,
                measurement_name,
                X_layer_name="data",
                extra_X_layer_names=["nonesuch"],
            )

        bdata = tiledbsoma.io.to_anndata(
            exp, measurement_name, X_layer_name="data", extra_X_layer_names=["data2"]
        )
        assert bdata.X is not None
        assert len(bdata.layers) == 1
        assert sorted(list(bdata.layers.keys())) == ["data2"]

        # extra_X_layer_names as list
        bdata = tiledbsoma.io.to_anndata(
            exp,
            measurement_name,
            X_layer_name="data",
            extra_X_layer_names=["data2", "data"],
        )
        assert bdata.X is not None
        assert len(bdata.layers) == 1
        assert sorted(list(bdata.layers.keys())) == ["data2"]

        # extra_X_layer_names as tuple
        bdata = tiledbsoma.io.to_anndata(
            exp,
            measurement_name,
            X_layer_name="data",
            extra_X_layer_names=("data2", "data", "data3"),
        )
        assert bdata.X is not None
        assert len(bdata.layers) == 2
        assert sorted(list(bdata.layers.keys())) == ["data2", "data3"]

        # measurement.X.keys() is a collections.abc.KeysView -- check that they can do this without
        # having to do list(measurement.X.keys())
        bdata = tiledbsoma.io.to_anndata(
            exp,
            measurement_name,
            X_layer_name="data",
            extra_X_layer_names=measurement.X.keys(),
        )
        assert bdata.X is not None
        assert len(bdata.layers) == 2
        assert sorted(list(bdata.layers.keys())) == ["data2", "data3"]


# fmt: off
@pytest.mark.parametrize("dtype", ["float64", "string"])          # new column dtype
@pytest.mark.parametrize("nans", ["all", "none", "some"])         # how many `nan`s in new column?
@pytest.mark.parametrize("new_obs_ids", ["all", "none", "half"])  # how many new obs IDs?
# fmt: on
def test_nan_append(conftest_pbmc_small, dtype, nans, new_obs_ids):
    """Test append-ingesting an AnnData object, including a new `obs` column with various properties:

    - {all,some,none} of its values are `nan`
    - dtype string or float
    - {all,none,half} of its obs IDs already present in the 1st AnnData.

    Prompted by observing a failure when all obs IDs are already present, and a new column has string dtype and contains
    at least one `nan`. See also:
    - https://github.com/single-cell-data/TileDB-SOMA/pull/2357
    - https://github.com/single-cell-data/TileDB-SOMA/pull/2364
    """
    conftest_pbmc_small.obsm = None
    conftest_pbmc_small.varm = None
    conftest_pbmc_small.obsp = None
    conftest_pbmc_small.varp = None
    conftest_pbmc_small.uns = dict()

    # Add empty column to obs
    obs = conftest_pbmc_small.obs
    if nans == "all":
        obs["batch_id"] = np.nan
    else:
        elem = "batch_id" if dtype == "string" else 1.23
        if nans == "some":
            obs["batch_id"] = np.nan
            obs.loc[obs.index.tolist()[0], "batch_id"] = elem
        else:
            obs["batch_id"] = elem

    obs["batch_id"] = obs["batch_id"].astype(dtype)

    # Create a copy of the anndata object
    adata2 = conftest_pbmc_small.copy()
    obs2 = adata2.obs
    if new_obs_ids == "all":
        obs2.index = obs2.index + "_2"
    elif new_obs_ids == "half":
        half = len(obs2) // 2
        obs2.index = obs2.index[:half].tolist() + (obs2.index[half:] + "_2").tolist()

    # Initial ingest
    SOMA_URI = tempfile.mkdtemp(prefix="soma-exp-")
    tiledbsoma.io.from_anndata(
        experiment_uri=SOMA_URI, anndata=conftest_pbmc_small, measurement_name="RNA"
    )

    # Register the second anndata object
    rd = tiledbsoma.io.register_anndatas(
        experiment_uri=SOMA_URI,
        adatas=[adata2],
        measurement_name="RNA",
        obs_field_name="obs_id",
        var_field_name="var_id",
    )

    nobs = rd.get_obs_shape()
    nvars = rd.get_var_shapes()
    tiledbsoma.io.resize_experiment(SOMA_URI, nobs=nobs, nvars=nvars)

    # Append the second anndata object
    tiledbsoma.io.from_anndata(
        experiment_uri=SOMA_URI,
        anndata=adata2,
        measurement_name="RNA",
        registration_mapping=rd,
    )


# This tests https://github.com/single-cell-data/TileDB-SOMA/pull/3354:
# * High-cardinality categorical-of-string gets converted to plain string, for
#   a significant performance benefit in subsequent processing
# * Two particular case needing testing:
#   o When a high-cardinality enum is appended to existing storage with categorical of string;
#   o When a low-cardinality enum is appended to existing storage with plain string.
def test_decat_append(tmp_path):

    # Prepare the AnnData inputs
    nobs_under = tiledbsoma.io.conversions.COLUMN_DECAT_THRESHOLD - 2
    nobs_over = tiledbsoma.io.conversions.COLUMN_DECAT_THRESHOLD + 2
    nvar = 100

    obs_ids_under = [f"under_{e:08}" for e in range(nobs_under)]
    obs_ids_over = [f"over_{e:08}" for e in range(nobs_over)]
    var_ids = [f"gene_{e:08}" for e in range(nvar)]

    string_enum_values_under = [f"enum_u_{e:06}" for e in range(nobs_under)]
    string_enum_values_over = [f"enum_o_{e:06}" for e in range(nobs_over)]
    float_enum_values_under = [1e6 + e for e in range(nobs_under)]
    float_enum_values_over = [2e6 + e for e in range(nobs_over)]
    bool_enum_values_under = random.choices([True, False], k=nobs_under)
    bool_enum_values_over = random.choices([True, False], k=nobs_over)

    obs_under = pd.DataFrame(
        data={
            "obs_id": np.asarray(obs_ids_under),
            "is_primary_data": np.asarray([True] * nobs_under),
            "string_enum": pd.Series(
                np.asarray(string_enum_values_under), dtype="category"
            ),
            "float_enum": pd.Series(
                np.asarray(float_enum_values_under), dtype="category"
            ),
            "bool_enum": pd.Series(
                np.asarray(bool_enum_values_under), dtype="category"
            ),
        }
    )
    obs_under.set_index("obs_id", inplace=True)

    obs_over = pd.DataFrame(
        data={
            "obs_id": np.asarray(obs_ids_over),
            "is_primary_data": np.asarray([True] * nobs_over),
            "string_enum": pd.Series(
                np.asarray(string_enum_values_over), dtype="category"
            ),
            "float_enum": pd.Series(
                np.asarray(float_enum_values_over), dtype="category"
            ),
            "bool_enum": pd.Series(np.asarray(bool_enum_values_over), dtype="category"),
        }
    )
    obs_over.set_index("obs_id", inplace=True)

    var = pd.DataFrame(
        data={
            "var_id": np.asarray(var_ids),
            "mybool": np.asarray([True] * nvar),
        }
    )
    var.set_index("var_id", inplace=True)

    X_under = scipy.sparse.random(
        nobs_under, nvar, density=0.1, dtype=np.float64
    ).tocsr()
    X_over = scipy.sparse.random(nobs_over, nvar, density=0.1, dtype=np.float64).tocsr()

    adata_under = anndata.AnnData(
        X=X_under, obs=obs_under, var=var, dtype=X_under.dtype
    )
    adata_over = anndata.AnnData(X=X_over, obs=obs_over, var=var, dtype=X_over.dtype)

    # Do the initial ingests from AnnData format to TileDB-SOMA format
    path_under = (tmp_path / "under").as_posix()
    path_over = (tmp_path / "over").as_posix()

    tiledbsoma.io.from_anndata(path_under, adata_under, "RNA")
    tiledbsoma.io.from_anndata(path_over, adata_over, "RNA")

    # Check that the low-cardinality categorical-of-string in the AnnData has
    # been ingested to TileDB-SOMA enum-of-string.
    with tiledbsoma.Experiment.open(path_under) as exp_under:
        assert pa.types.is_dictionary(exp_under.obs.schema.field("string_enum").type)
        assert pa.types.is_dictionary(exp_under.obs.schema.field("float_enum").type)
        assert pa.types.is_dictionary(exp_under.obs.schema.field("bool_enum").type)
        obs_table = exp_under.obs.read().concat()
        assert obs_table.column("string_enum").to_pylist() == string_enum_values_under
        assert obs_table.column("float_enum").to_pylist() == float_enum_values_under
        assert obs_table.column("bool_enum").to_pylist() == bool_enum_values_under

    # Check that the high-cardinality categorical-of-string in the AnnData has
    # been ingested to TileDB-SOMA plain string.
    with tiledbsoma.Experiment.open(path_over) as exp_over:
        assert not pa.types.is_dictionary(exp_over.obs.schema.field("string_enum").type)
        assert not pa.types.is_dictionary(exp_over.obs.schema.field("float_enum").type)
        assert pa.types.is_dictionary(exp_over.obs.schema.field("bool_enum").type)
        obs_table = exp_over.obs.read().concat()
        assert obs_table.column("string_enum").to_pylist() == string_enum_values_over
        assert obs_table.column("float_enum").to_pylist() == float_enum_values_over
        assert obs_table.column("bool_enum").to_pylist() == bool_enum_values_over

    # Append over-the-threshold AnnData to under-the-threshold TileDB-SOMA
    # storage, and vice versa.
    rd_under_over = tiledbsoma.io.register_anndatas(
        experiment_uri=path_under,
        adatas=[adata_over],
        measurement_name="RNA",
        obs_field_name="obs_id",
        var_field_name="var_id",
    )

    rd_over_under = tiledbsoma.io.register_anndatas(
        experiment_uri=path_over,
        adatas=[adata_under],
        measurement_name="RNA",
        obs_field_name="obs_id",
        var_field_name="var_id",
    )

    tiledbsoma.io.resize_experiment(
        path_under,
        nobs=rd_under_over.get_obs_shape(),
        nvars=rd_under_over.get_var_shapes(),
    )
    tiledbsoma.io.resize_experiment(
        path_over,
        nobs=rd_over_under.get_obs_shape(),
        nvars=rd_over_under.get_var_shapes(),
    )

    tiledbsoma.io.from_anndata(
        path_under, adata_over, "RNA", registration_mapping=rd_under_over
    )
    tiledbsoma.io.from_anndata(
        path_over, adata_under, "RNA", registration_mapping=rd_over_under
    )

    # Check that the appends happened successfully
    with tiledbsoma.Experiment.open(path_under) as exp_under:
        assert pa.types.is_dictionary(exp_under.obs.schema.field("string_enum").type)
        assert pa.types.is_dictionary(exp_under.obs.schema.field("float_enum").type)
        assert pa.types.is_dictionary(exp_under.obs.schema.field("bool_enum").type)
        obs_table = exp_under.obs.read().concat()
        assert (
            obs_table.column("string_enum").to_pylist()
            == string_enum_values_under + string_enum_values_over
        )
        assert (
            obs_table.column("float_enum").to_pylist()
            == float_enum_values_under + float_enum_values_over
        )
        assert (
            obs_table.column("bool_enum").to_pylist()
            == bool_enum_values_under + bool_enum_values_over
        )

    with tiledbsoma.Experiment.open(path_over) as exp_over:
        assert not pa.types.is_dictionary(exp_over.obs.schema.field("string_enum").type)
        assert not pa.types.is_dictionary(exp_over.obs.schema.field("float_enum").type)
        assert pa.types.is_dictionary(exp_over.obs.schema.field("bool_enum").type)
        obs_table = exp_over.obs.read().concat()
        assert (
            obs_table.column("string_enum").to_pylist()
            == string_enum_values_over + string_enum_values_under
        )
        assert (
            obs_table.column("float_enum").to_pylist()
            == float_enum_values_over + float_enum_values_under
        )
        assert (
            obs_table.column("bool_enum").to_pylist()
            == bool_enum_values_over + bool_enum_values_under
        )


def test_from_h5ad_bad_uri():
    with pytest.raises(tiledbsoma.SOMAError, match="URI /nonesuch is not a valid URI"):
        next(tiledbsoma.io._util.read_h5ad("/nonesuch").gen)
