import pathlib
import shutil

import pytest

import tiledbsoma as soma

from ._util import ROOT_DATA_DIR
from .test_experiment_query_spatial import soma_spatial_experiment  # noqa: F401


@pytest.mark.parametrize(
    "experiment_path",
    [
        ROOT_DATA_DIR / "soma-experiment-versions-2025-04-04" / ver / "pbmc3k_processed"
        for ver in ("1.12.3", "1.14.5", "1.15.0", "1.15.7", "1.7.3")
    ],
)
@pytest.mark.parametrize(
    "coords",
    [
        (slice(None),),  # all
        (),  # none
        (slice(10, 100),),
        (list(range(0, 1000, 7)),),
    ],
)
@pytest.mark.parametrize(
    "value_filter",
    [
        None,  # no filter
        "n_genes > 500",
        "louvain == 'B cells' and n_genes < 500",
    ],
)
def test_experiment_obs_axis_delete_from_pbmc3k(
    tmp_path, soma_tiledb_context, experiment_path, coords, value_filter
) -> None:
    """Test on PBMC3K dataset of various vintages."""

    # Make a copy of the Experiment as to not write over the original
    exp_path = pathlib.PosixPath(experiment_path)
    uri = (tmp_path / exp_path.name).as_posix()
    shutil.copytree(exp_path, uri)

    # Grab the list of joinids we expect to be deleted
    with soma.Experiment.open(uri, mode="r", context=soma_tiledb_context) as exp:
        joinids = (
            exp.obs.read(
                coords=coords,
                value_filter=value_filter,
                column_names=["soma_joinid"],
            )
            .concat()
            .column("soma_joinid")
            .combine_chunks()
        )
    assert len(joinids)  # just in case we mess up the test params

    # Do the delete
    with soma.Experiment.open(uri, mode="d", context=soma_tiledb_context) as exp:
        exp.obs_axis_delete(coords, value_filter=value_filter)

    # Verify the result
    with soma.Experiment.open(uri, mode="r", context=soma_tiledb_context) as exp:
        assert len(exp.obs.read(coords=(joinids,)).concat()) == 0
        for ms in exp.ms.values():
            for arr in ms.X.values():
                assert len(arr.read(coords=(joinids,)).tables().concat()) == 0
            if "obsm" in ms:
                for arr in ms.obsm.values():
                    assert len(arr.read(coords=(joinids,)).tables().concat()) == 0
            if "obsp" in ms:
                for arr in ms.obsp.values():
                    assert len(arr.read(coords=(joinids, slice(None))).tables().concat()) == 0
                    assert len(arr.read(coords=(slice(None), joinids)).tables().concat()) == 0

        # we don't expect these in this dataset
        assert "obs_spatial_presence" not in exp
        assert "spatial" not in exp


# var_id, n_cells
# len  1837


@pytest.mark.parametrize(
    "experiment_path",
    [
        ROOT_DATA_DIR / "soma-experiment-versions-2025-04-04" / ver / "pbmc3k_processed"
        for ver in ("1.12.3", "1.14.5", "1.15.0", "1.15.7", "1.7.3")
    ],
)
@pytest.mark.parametrize(
    "coords",
    [
        (slice(None),),  # all
        (),  # none
        (slice(10, 100),),
        (list(range(0, 1000, 7)),),
    ],
)
@pytest.mark.parametrize(
    "value_filter",
    [
        None,  # no filter
        "n_cells > 500",
        "var_id in ['PRDX1', 'TMED5', 'CDA', 'C1QA', 'C1QC', 'C1QB', 'ZNF436']",
    ],
)
def test_experiment_var_axis_delete_from_pbmc3k(
    tmp_path, soma_tiledb_context, experiment_path, coords, value_filter
) -> None:
    """Test on PBMC3K dataset of various vintages."""

    # Make a copy of the Experiment as to not write over the original
    exp_path = pathlib.PosixPath(experiment_path)
    uri = (tmp_path / exp_path.name).as_posix()
    shutil.copytree(exp_path, uri)

    measurement_name = "RNA"

    # Grab the list of joinids we expect to be deleted
    with soma.Experiment.open(uri, mode="r", context=soma_tiledb_context) as exp:
        joinids = (
            exp.ms[measurement_name]
            .var.read(
                coords=coords,
                value_filter=value_filter,
                column_names=["soma_joinid"],
            )
            .concat()
            .column("soma_joinid")
            .combine_chunks()
        )
    assert len(joinids)  # just in case we mess up the test params

    # Do the delete
    with soma.Experiment.open(uri, mode="d", context=soma_tiledb_context) as exp:
        exp.var_axis_delete(measurement_name, coords, value_filter=value_filter)

    # Verify the result
    with soma.Experiment.open(uri, mode="r", context=soma_tiledb_context) as exp:
        assert len(exp.ms[measurement_name].var.read(coords=(joinids,)).concat()) == 0
        for arr in exp.ms[measurement_name].X.values():
            assert len(arr.read(coords=(slice(None), joinids)).tables().concat()) == 0
        if "varm" in exp.ms[measurement_name]:
            for arr in exp.ms[measurement_name].varm.values():
                assert len(arr.read(coords=(joinids,)).tables().concat()) == 0
        if "varp" in exp.ms[measurement_name]:
            for arr in exp.ms[measurement_name].varp.values():
                assert len(arr.read(coords=(joinids, slice(None))).tables().concat()) == 0
                assert len(arr.read(coords=(slice(None), joinids)).tables().concat()) == 0

        # we don't expect these in this dataset
        assert "var_spatial_presence" not in exp.ms[measurement_name]
        assert "spatial" not in exp

    # ALSO verify we did not touch the other measurements
    with soma.open(experiment_path.as_posix(), context=soma_tiledb_context) as orig_exp:
        for ms_name in exp.ms:
            if ms_name == measurement_name:
                continue
            assert exp.ms[ms_name].var.count == orig_exp.ms[ms_name].var.count
            for arr_name in exp.ms[ms_name].X:
                assert exp.ms[ms_name].X[arr_name].nnz == orig_exp.ms[ms_name].X[arr_name].nnz


@pytest.mark.spatial
def test_experiment_obs_axis_delete_spatial(soma_spatial_experiment, soma_tiledb_context) -> None:  # noqa: F811
    with soma.open(soma_spatial_experiment.uri, mode="d", context=soma_tiledb_context) as exp:
        exp.obs_axis_delete(value_filter="soma_joinid >= 33")

    exp.reopen(mode="r")
    assert (exp.obs.read().concat().to_pandas()["soma_joinid"] < 33).all()
    assert (exp.obs_spatial_presence.read().concat().to_pandas()["soma_joinid"] < 33).all()
    for ms in exp.ms.values():
        for arr in ms.X.values():
            assert (arr.read().tables().concat().to_pandas()["soma_dim_0"] < 33).all()
    for sc in exp.spatial.values():
        for arr in sc.obsl.values():
            assert (arr.read().concat().to_pandas()["soma_joinid"] < 33).all()


"""Other test TODO:
1. make an experiment with all slots, and make sure we delete from them (mock delete_cells)
2. check that we catch dense arrays and error out
3. check that we ignore weird stuff
"""
