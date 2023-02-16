import tempfile
from pathlib import Path

import anndata
import pytest
import tiledb

import tiledbsoma
import tiledbsoma.io
from tiledbsoma import _constants, _factory

HERE = Path(__file__).parent


@pytest.fixture
def h5ad_file(request):
    # pbmc-small is faster for automated unit-test / CI runs.
    # input_path = HERE.parent / "testdata/pbmc3k_processed.h5ad"
    input_path = HERE.parent / "testdata/pbmc-small.h5ad"
    return input_path


@pytest.fixture
def h5ad_file_extended(request):
    # This has more component arrays in it
    input_path = HERE.parent / "testdata/pbmc3k_processed.h5ad"
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
def test_import_anndata(adata, ingest_modes):

    have_ingested = False

    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name
    orig = adata
    metakey = _constants.SOMA_OBJECT_TYPE_METADATA_KEY  # keystroke-saver
    all2d = (slice(None), slice(None))  # keystroke-saver

    for ingest_mode in ingest_modes:

        write_exp = tiledbsoma.io.from_anndata(
            output_path, orig, "RNA", ingest_mode=ingest_mode
        )
        if ingest_mode != "schema_only":
            have_ingested = True
        write_exp.close()
        # del write_exp

        exp = _factory.open(output_path)
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

        # Check X/data (dense)
        assert exp.ms["RNA"].X["data"].metadata.get(metakey) == "SOMADenseNDArray"
        if have_ingested:
            matrix = exp.ms["RNA"].X["data"].read(coords=all2d)
            assert matrix.shape == orig.X.shape
        else:
            with pytest.raises(ValueError):
                exp.ms["RNA"].X["data"].read(coords=all2d)

        # Check raw/X/data (sparse)
        assert exp.ms["raw"].X["data"].metadata.get(metakey) == "SOMASparseNDArray"
        if have_ingested:
            matrix = exp.ms["raw"].X["data"].read(coords=all2d).coos().concat()
            assert matrix.shape == orig.raw.X.shape

        # Check some annotation matrices

        obsm = exp.ms["RNA"].obsm
        assert sorted(obsm.keys()) == sorted(orig.obsm.keys())
        for key in list(orig.obsm.keys()):
            assert obsm[key].metadata.get(metakey) == "SOMADenseNDArray"
            if have_ingested:
                matrix = obsm[key].read(coords=all2d)
                assert matrix.shape == orig.obsm[key].shape
            else:
                with pytest.raises(ValueError):
                    matrix = obsm[key].read(coords=all2d)

        varm = exp.ms["RNA"].varm
        assert sorted(varm.keys()) == sorted(orig.varm.keys())
        for key in list(orig.varm.keys()):
            assert varm[key].metadata.get(metakey) == "SOMADenseNDArray"
            if have_ingested:
                matrix = varm[key].read(coords=all2d)
                assert matrix.shape == orig.varm[key].shape
            else:
                with pytest.raises(ValueError):
                    matrix = varm[key].read(coords=all2d)

        obsp = exp.ms["RNA"].obsp
        assert sorted(obsp.keys()) == sorted(orig.obsp.keys())
        for key in list(orig.obsp.keys()):
            assert obsp[key].metadata.get(metakey) == "SOMASparseNDArray"
            if have_ingested:
                matrix = obsp[key].read(coords=all2d).coos().concat()
                assert matrix.shape == orig.obsp[key].shape

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
        output_path1, resume_mode_h5ad_file, "RNA", ingest_mode="write"
    ).close()

    tempdir2 = tempfile.TemporaryDirectory()
    output_path2 = tempdir2.name
    start_write = tiledbsoma.io.from_h5ad(
        output_path2, resume_mode_h5ad_file, "RNA", ingest_mode="write"
    )
    start_write.close()
    tiledbsoma.io.from_h5ad(
        output_path2, resume_mode_h5ad_file, "RNA", ingest_mode="resume"
    ).close()

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
        h5ad_file_extended,
        measurement_name="RNA",
        use_relative_uri=use_relative_uri,
    ).close()

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


def test_add_matrix_to_collection(adata):
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name

    exp = tiledbsoma.io.from_anndata(output_path, adata, measurement_name="RNA")
    exp.close()
    with _factory.open(output_path) as exp_r:
        assert list(exp_r.ms["RNA"].X.keys()) == ["data"]

    with _factory.open(output_path, "w") as exp:
        tiledbsoma.io.add_X_layer(exp, "RNA", "data2", adata.X)
    with _factory.open(output_path) as exp_r:
        assert sorted(list(exp_r.ms["RNA"].X.keys())) == ["data", "data2"]

    with _factory.open(output_path, "w") as exp:
        with pytest.raises(KeyError):
            tiledbsoma.io.add_X_layer(exp, "nonesuch", "data3", adata.X)

    with _factory.open(output_path) as exp_r:
        assert sorted(list(exp_r.ms["RNA"].obsm.keys())) == ["X_pca", "X_tsne"]
    with _factory.open(output_path, "w") as exp:
        tiledbsoma.io.add_matrix_to_collection(
            exp, "RNA", "obsm", "X_pcb", adata.obsm["X_pca"]
        )
    with _factory.open(output_path) as exp_r:
        assert sorted(list(exp_r.ms["RNA"].obsm.keys())) == ["X_pca", "X_pcb", "X_tsne"]

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


def test_export_anndata(adata):
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name

    exp = tiledbsoma.io.from_anndata(output_path, adata, measurement_name="RNA")
    exp.close()

    with _factory.open(output_path) as exp:
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
