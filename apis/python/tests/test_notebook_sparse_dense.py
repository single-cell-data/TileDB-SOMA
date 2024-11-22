import tarfile

import pytest

import tiledbsoma
import tiledbsoma.io

from tests._util import PY_ROOT


@pytest.mark.parametrize("name", ["sparse", "dense"])
def test_notebook_path_dense(tmp_path, name):
    tgz_path = PY_ROOT / f"notebooks/data/pbmc3k-{name}.tgz"
    uri = tmp_path.as_posix()

    with tarfile.open(tgz_path) as handle:
        handle.extractall(uri)

    with tiledbsoma.Experiment.open(uri) as exp:
        assert len(exp.obs.read().concat()) == 2638
