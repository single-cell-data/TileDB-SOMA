import pytest

import tiledbsoma
import tiledbsoma.io

from tests._util import PY_ROOT


@pytest.mark.parametrize("name", ["sparse", "dense"])
def test_notebook_path_dense(name):
    path = PY_ROOT / f"notebooks/data/{name}/pbmc3k"

    with tiledbsoma.Experiment.open(path.as_posix()) as exp:
        assert len(exp.obs.read().concat()) == 2638
