import os

import pytest
from typeguard import install_import_hook

# Must be called before tiledbsoma imports
install_import_hook("tiledbsoma")

import tiledbsoma  # noqa: E402
import tiledb.cloud  # noqa: E402


@pytest.fixture
def conftest_token():
    env_name = "TILEDB_REST_UNITTEST_TOKEN"
    token = os.getenv(env_name)
    if token is None:
        raise Exception(f'Environment variable "{env_name}" is not set')
    return token


@pytest.fixture
def conftest_tiledb_cloud_login(conftest_token):
    print("conftest_tiledb_cloud_login")
    tiledb.cloud.login(token=conftest_token)
    return


@pytest.fixture
def conftest_user_profile(conftest_tiledb_cloud_login):
    return tiledb.cloud.user_profile()


@pytest.fixture
def conftest_namespace(conftest_user_profile):
    return conftest_user_profile.username


@pytest.fixture
def conftest_default_s3_path(conftest_user_profile):
    return conftest_user_profile.default_s3_path


@pytest.fixture
def conftest_context(conftest_token, conftest_namespace):
    return tiledbsoma.SOMAContext(
        config={
            "rest.token": conftest_token,
            "rest.payer_namespace": conftest_namespace,
        },
    )
