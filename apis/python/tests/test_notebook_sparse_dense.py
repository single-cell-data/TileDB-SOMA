from pathlib import Path

import pytest

import tiledbsoma
import tiledbsoma.io

HERE = Path(__file__).parent


@pytest.mark.parametrize("name", ["sparse", "dense"])
def test_notebook_path_dense(name):
    path = HERE.parent / f"notebooks/data/{name}/pbmc3k"

    with tiledbsoma.Experiment.open(path.as_posix()) as exp:
        assert len(exp.obs.read().concat()) == 2638
