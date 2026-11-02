# Tests that run a UDF on TileDB Cloud.
#
# UDFs are pickled locally and unpickled server-side, so the local Python version must match
# the TileDB Cloud UDF runtime. That can be older than the oldest Python tiledbsoma supports,
# so these tests deliberately live outside remote_tests/ (whose conftest imports tiledbsoma)
# and must not import tiledbsoma locally. Inside the UDF, `import tiledbsoma` resolves to the
# version installed on TileDB Cloud.

from __future__ import annotations

import os
import sys

import pytest

import tiledb.cloud

if os.getenv("TILEDB_REST_UNITTEST_TOKEN") is None:
    pytest.skip(reason="$TILEDB_REST_UNITTEST_TOKEN is not set", allow_module_level=True)


@pytest.fixture
def tiledb_cloud_login():
    tiledb.cloud.login(token=os.environ["TILEDB_REST_UNITTEST_TOKEN"])


@pytest.mark.skipif(
    (sys.version_info.major, sys.version_info.minor) != (3, 9),
    reason="As of 2025-02-05 UDFs require Python 3.9",
)
def test_remote_version(tiledb_cloud_login):
    def remote_version():
        import tiledbsoma

        return {"tiledbsoma": tiledbsoma.__version__}

    output = tiledb.cloud.udf.exec(remote_version)
    assert "tiledbsoma" in output
    assert output["tiledbsoma"].startswith("2.")
