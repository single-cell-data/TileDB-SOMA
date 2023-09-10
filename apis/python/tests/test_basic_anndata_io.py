import pathlib
import tempfile
from pathlib import Path
from typing import Optional

import anndata
import numpy as np
import pandas as pd
import pytest
import somacore
import tiledb

import tiledbsoma
import tiledbsoma.io
from tiledbsoma import _constants, _factory

HERE = Path(__file__).parent


@pytest.fixture
def h5ad_file(request):
    # pbmc-small is faster for automated unit-test / CI runs.
    input_path = HERE.parent / "testdata/pbmc-small.h5ad"
    # input_path = HERE.parent / "testdata/pbmc3k_processed.h5ad"
    return input_path


@pytest.fixture
def h5ad_file_extended(request):
    # This has more component arrays in it
    input_path = HERE.parent / "testdata/pbmc3k_processed.h5ad"
    return input_path


@pytest.fixture
def h5ad_file_with_obsm_holes(request):
    # This has zeroes in an obsm matrix so nnz is not num_rows * num_cols
    input_path = HERE.parent / "testdata/pbmc3k-with-obsm-zero.h5ad"
    return input_path


@pytest.fixture
def h5ad_file_uns_string_array(request):
    # This has uns["louvain_colors"] with dtype.char == "U"
    input_path = HERE.parent / "testdata/pbmc3k.h5ad"
    return input_path


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
    input_path = HERE.parent / "testdata/categorical_int_nan.h5ad"
    return input_path


@pytest.fixture
def adata(h5ad_file):
    return anndata.read_h5ad(h5ad_file)


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
def test_import_anndata(adata, ingest_modes, X_kind):
    have_ingested = False

    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name
    orig = adata
    metakey = _constants.SOMA_OBJECT_TYPE_METADATA_KEY  # keystroke-saver
    all2d = (slice(None), slice(None))  # keystroke-saver

    adata.layers["plus1"] = adata.X + 1

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
                # Only sparse arrays have used_shape
                if ingest_mode != "schema_only":
                    # E.g. AnnData X shape is (80,20). We expect used_shape ((0,79),(0,19)).
                    assert X.used_shape() == tuple([(0, e - 1) for e in orig.shape])
            else:
                assert X.metadata.get(metakey) == "SOMADenseNDArray"
                if have_ingested:
                    matrix = X.read(coords=all2d)
                    assert matrix.size == orig.X.size
                else:
                    with pytest.raises(ValueError):
                        X.read(coords=all2d)

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


def _get_fragment_count(array_uri):
    return len(tiledb.fragment.FragmentInfoList(array_uri=array_uri))


@pytest.mark.parametrize(
    "resume_mode_h5ad_file",
    [
        HERE.parent / "testdata/pbmc-small-x-dense.h5ad",
        HERE.parent / "testdata/pbmc-small-x-csr.h5ad",
        HERE.parent / "testdata/pbmc-small-x-csc.h5ad",
    ],
)
def test_resume_mode(adata, resume_mode_h5ad_file):
    """
    Makes sure resume-mode ingest after successful ingest of the same input data does not write
    anything new
    """

    tempdir1 = tempfile.TemporaryDirectory()
    output_path1 = tempdir1.name
    tiledbsoma.io.from_h5ad(
        output_path1, resume_mode_h5ad_file.as_posix(), "RNA", ingest_mode="write"
    )

    tempdir2 = tempfile.TemporaryDirectory()
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
def test_ingest_relative(h5ad_file_extended, use_relative_uri):
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name

    tiledbsoma.io.from_h5ad(
        output_path,
        h5ad_file_extended.as_posix(),
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

    exp = tiledbsoma.open(output_path)
    with tiledb.Group(exp.uri) as G:
        assert G.is_relative("obs") == expected_relative
        assert G.is_relative("ms") == expected_relative

    with tiledb.Group(exp.ms.uri) as G:
        assert G.is_relative("RNA") == expected_relative
    with tiledb.Group(exp.ms["RNA"].uri) as G:
        assert G.is_relative("var") == expected_relative
        assert G.is_relative("X") == expected_relative
    with tiledb.Group(exp.ms["RNA"].X.uri) as G:
        assert G.is_relative("data") == expected_relative

    for collection_name in ["obsm", "obsp", "varm"]:  # h5ad_file_extended has no varp
        with tiledb.Group(exp.ms["RNA"][collection_name].uri) as G:
            for member in G:
                assert G.is_relative(member.name) == expected_relative

    with tiledb.Group(exp.ms.uri) as G:
        assert G.is_relative("raw") == expected_relative
    with tiledb.Group(exp.ms["raw"].uri) as G:
        assert G.is_relative("var") == expected_relative
        assert G.is_relative("X") == expected_relative
    with tiledb.Group(exp.ms["raw"].X.uri) as G:
        assert G.is_relative("data") == expected_relative

    exp.close()


def test_ingest_uns(tmp_path: pathlib.Path, h5ad_file_extended):
    tmp_uri = tmp_path.as_uri()
    original = anndata.read(h5ad_file_extended)
    uri = tiledbsoma.io.from_anndata(tmp_uri, original, measurement_name="hello")

    with tiledbsoma.Experiment.open(uri) as exp:
        uns = exp.ms["hello"]["uns"]
        assert isinstance(uns, tiledbsoma.Collection)
        assert uns.metadata["soma_tiledbsoma_type"] == "uns"
        assert set(uns) == {
            "draw_graph",
            "louvain",
            "louvain_colors",
            "neighbors",
            "pca",
            "rank_genes_groups",
        }
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
        assert np.array_equal(got_pca_variance, original.uns["pca"]["variance"])


def test_ingest_uns_string_array(h5ad_file_uns_string_array):
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name

    tiledbsoma.io.from_h5ad(
        output_path,
        h5ad_file_uns_string_array.as_posix(),
        measurement_name="RNA",
    )

    with tiledbsoma.Experiment.open(output_path) as exp:
        with tiledbsoma.DataFrame.open(
            exp.ms["RNA"]["uns"]["louvain_colors"].uri
        ) as df:
            contents = df.read().concat()["values"]
            assert len(contents) == 8
            assert contents[0].as_py() == "#1f77b4"


def test_add_matrix_to_collection(adata):
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name

    uri = tiledbsoma.io.from_anndata(output_path, adata, measurement_name="RNA")
    exp = tiledbsoma.Experiment.open(uri)
    with _factory.open(output_path) as exp_r:
        assert list(exp_r.ms["RNA"].X.keys()) == ["data"]
        with pytest.raises(tiledbsoma.SOMAError):
            tiledbsoma.io.add_X_layer(exp, "RNA", "data2", adata.X)  # not open for read
    with _factory.open(output_path, "w") as exp:
        tiledbsoma.io.add_X_layer(exp, "RNA", "data2", adata.X)
    with pytest.raises(tiledbsoma.SOMAError):
        tiledbsoma.io.add_X_layer(exp, "RNA", "data3", adata.X)  # closed
    with _factory.open(output_path) as exp_r:
        assert sorted(list(exp_r.ms["RNA"].X.keys())) == ["data", "data2"]

    with _factory.open(output_path, "w") as exp:
        with pytest.raises(KeyError):
            tiledbsoma.io.add_X_layer(exp, "nonesuch", "data3", adata.X)

    with _factory.open(output_path) as exp_r:
        assert sorted(list(exp_r.ms["RNA"].obsm.keys())) == sorted(
            list(adata.obsm.keys())
        )
    with _factory.open(output_path, "w") as exp:
        tiledbsoma.io.add_matrix_to_collection(
            exp, "RNA", "obsm", "X_pcb", adata.obsm["X_pca"]
        )
    with _factory.open(output_path) as exp_r:
        assert sorted(list(exp_r.ms["RNA"].obsm.keys())) == sorted(
            list(adata.obsm.keys()) + ["X_pcb"]
        )

    with _factory.open(output_path, "w") as exp:
        with pytest.raises(KeyError):
            tiledbsoma.io.add_matrix_to_collection(
                exp, "nonesuch", "obsm", "X_pcc", adata.obsm["X_pca"]
            )

        tiledbsoma.io.add_matrix_to_collection(
            exp, "RNA", "newthing", "X_pcd", adata.obsm["X_pca"]
        )
    with _factory.open(output_path) as exp_r:
        assert sorted(list(exp_r.ms["RNA"]["newthing"].keys())) == ["X_pcd"]


# This tests tiledbsoma.io._create_or_open_coll -> tiledbsoma.io._create_or_open_collection
# compatibility. As a rule we do not offer any backward-compatibility assurances for non-public
# modules or methods -- those whose names start with an underscore.  For this single case we are
# making an exception. For future code-imitation purposes, please be aware this is a pattern to be
# avoided in the future, not imitated.
def test_add_matrix_to_collection_1_2_7(adata):
    def add_X_layer(
        exp: tiledbsoma.Experiment,
        measurement_name: str,
        X_layer_name: str,
        X_layer_data,  # E.g. a scipy.csr_matrix from scanpy analysis
        ingest_mode: str = "write",
        use_relative_uri: Optional[bool] = None,
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
        use_relative_uri: Optional[bool] = None,
        context: Optional[tiledbsoma.SOMATileDBContext] = None,
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

    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name

    uri = tiledbsoma.io.from_anndata(output_path, adata, measurement_name="RNA")
    exp = tiledbsoma.Experiment.open(uri)
    with _factory.open(output_path) as exp_r:
        assert list(exp_r.ms["RNA"].X.keys()) == ["data"]
        with pytest.raises(tiledbsoma.SOMAError):
            add_X_layer(exp, "RNA", "data2", adata.X)  # not open for read
    with _factory.open(output_path, "w") as exp:
        add_X_layer(exp, "RNA", "data2", adata.X)
    with pytest.raises(tiledbsoma.SOMAError):
        add_X_layer(exp, "RNA", "data3", adata.X)  # closed
    with _factory.open(output_path) as exp_r:
        assert sorted(list(exp_r.ms["RNA"].X.keys())) == ["data", "data2"]

    with _factory.open(output_path, "w") as exp:
        with pytest.raises(KeyError):
            add_X_layer(exp, "nonesuch", "data3", adata.X)

    with _factory.open(output_path) as exp_r:
        assert sorted(list(exp_r.ms["RNA"].obsm.keys())) == sorted(
            list(adata.obsm.keys())
        )

    with _factory.open(output_path, "w") as exp:
        add_matrix_to_collection(exp, "RNA", "obsm", "X_pcb", adata.obsm["X_pca"])
    with _factory.open(output_path) as exp_r:
        assert sorted(list(exp_r.ms["RNA"].obsm.keys())) == sorted(
            list(adata.obsm.keys()) + ["X_pcb"]
        )

    with _factory.open(output_path, "w") as exp:
        # It's nonsense biologically to add this to varp, but as a fake-test unit-test case, we can
        # use varp to test adding to a not-yet-existing collection.
        add_matrix_to_collection(exp, "RNA", "varp", "X_pcb", adata.obsm["X_pca"])
    with _factory.open(output_path) as exp_r:
        assert sorted(list(exp_r.ms["RNA"].varp.keys())) == sorted(
            list(adata.varp.keys()) + ["X_pcb"]
        )

    with _factory.open(output_path, "w") as exp:
        with pytest.raises(KeyError):
            add_matrix_to_collection(
                exp, "nonesuch", "obsm", "X_pcc", adata.obsm["X_pca"]
            )

        add_matrix_to_collection(exp, "RNA", "newthing", "X_pcd", adata.obsm["X_pca"])
    with _factory.open(output_path) as exp_r:
        assert sorted(list(exp_r.ms["RNA"]["newthing"].keys())) == ["X_pcd"]


def test_export_anndata(adata):
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name

    tiledbsoma.io.from_anndata(output_path, adata, measurement_name="RNA")

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

    assert readback.obs.shape == adata.obs.shape
    assert readback.var.shape == adata.var.shape
    assert readback.X.shape == adata.X.shape

    for key in adata.obsm.keys():
        assert readback.obsm[key].shape == adata.obsm[key].shape
    for key in adata.varm.keys():
        assert readback.varm[key].shape == adata.varm[key].shape
    for key in adata.obsp.keys():
        assert readback.obsp[key].shape == adata.obsp[key].shape
    for key in adata.varp.keys():
        assert readback.varp[key].shape == adata.varp[key].shape


def test_null_obs(adata, tmp_path: Path):
    output_path = tmp_path.as_uri()
    seed = 42
    #   Create column of all null values
    adata.obs["empty_all"] = pd.Categorical(
        [np.NaN] * adata.n_obs, dtype=pd.CategoricalDtype(categories=[], ordered=False)
    )
    #   Create column of partially-null values
    rng = np.random.RandomState(seed)
    adata.obs["empty_partial"] = rng.choice((np.NaN, 1.0), adata.n_obs, True)
    uri = tiledbsoma.io.from_anndata(
        output_path, adata, "RNA", ingest_mode="write", X_kind=tiledbsoma.SparseNDArray
    )
    exp = tiledbsoma.Experiment.open(uri)
    with tiledb.open(exp.obs.uri, "r") as obs:
        #   Explicitly check columns created above
        assert obs.attr("empty_all").isnullable
        assert obs.attr("empty_partial").isnullable
        #   For every column in the data frame
        #   ensure that `isnullable` reflects the null-ness
        #   of the Pandas data frame
        for k in adata.obs:
            assert obs.attr(k).isnullable == adata.obs[k].isnull().any()


# There exist in the wild AnnData files with categorical-int columns where the "not in the category"
# is indicated by the presence of floating-point math.NaN in cells. Here we test that we can ingest
# this.
def test_obs_with_categorical_int_nan_enumeration(
    tmp_path, h5ad_file_categorical_int_nan
):
    output_path = tmp_path.as_uri()

    tiledbsoma.io.from_h5ad(
        output_path, h5ad_file_categorical_int_nan, measurement_name="RNA"
    )

def test_export_obsm_with_holes(h5ad_file_with_obsm_holes, tmp_path):
    adata = anndata.read_h5ad(h5ad_file_with_obsm_holes.as_posix())
    assert 1 == 1

    # This data file is prepared such that obsm["X_pca"] has shape (2638, 50)
    # but its [0][0] element is a 0, so when it's stored as sparse, its nnz
    # is not 2638*50=131900.
    ado = adata.obsm["X_pca"]
    assert ado.shape == (2638, 50)

    output_path = tmp_path.as_posix()
    tiledbsoma.io.from_anndata(output_path, adata, "RNA")

    exp = tiledbsoma.Experiment.open(output_path)

    # Verify the bounding box on the SOMA SparseNDArray
    with tiledb.open(exp.ms["RNA"].obsm["X_pca"].uri) as so:
        assert so.meta["soma_dim_0_domain_lower"] == 0
        assert so.meta["soma_dim_0_domain_upper"] == 2637
        assert so.meta["soma_dim_1_domain_lower"] == 0
        assert so.meta["soma_dim_1_domain_upper"] == 49

    # With the bounding box present, all is well for outgest to AnnData format.
    try1 = tiledbsoma.io.to_anndata(exp, "RNA")
    assert try1.obsm["X_pca"].shape == (2638, 50)

    # Now remove the bounding box to simulate reading older data that lacks a bounding box.
    with tiledb.open(exp.ms["RNA"].obsm["X_pca"].uri, "w") as so:
        del so.meta["soma_dim_0_domain_lower"]
        del so.meta["soma_dim_0_domain_upper"]
        del so.meta["soma_dim_1_domain_lower"]
        del so.meta["soma_dim_1_domain_upper"]

    # Re-open to simulate opening afresh a bounding-box-free array.
    exp = tiledbsoma.Experiment.open(output_path)

    with tiledb.open(exp.ms["RNA"].obsm["X_pca"].uri) as so:
        with pytest.raises(KeyError):
            so.meta["soma_dim_0_domain_lower"]
        with pytest.raises(KeyError):
            so.meta["soma_dim_0_domain_upper"]
        with pytest.raises(KeyError):
            so.meta["soma_dim_1_domain_lower"]
        with pytest.raises(KeyError):
            so.meta["soma_dim_1_domain_upper"]
        assert so.meta["soma_object_type"] == "SOMASparseNDArray"

    # Now try the remaining options for outgest.
    with pytest.raises(tiledbsoma.SOMAError):
        tiledbsoma.io.to_anndata(exp, "RNA")

    try3 = tiledbsoma.io.to_anndata(
        exp, "RNA", obsm_varm_width_hints={"obsm": {"X_pca": 50}}
    )
    assert try3.obsm["X_pca"].shape == (2638, 50)