import datetime
import os
import pathlib
import shutil

import tiledb.cloud

# For cloud:
# * Create with timestamp
# * Delete on teardown
# For local:
# * Create without timestamp
#   o Only remove the URI from a _previous_ run (if any)
# * Do not delete on teardown -- so developers can look at the data


def util_make_uri(
    dirname: str,
    basename: str,
    namespace: str,
    default_s3_path: str,
) -> tuple[str, str]:
    if os.getenv("TILEDB_SOMA_CLOUD_TEST_LOCAL_PATHS") is None:
        # The default_s3_path contains the "s3://..." prefix and a trailing slash.
        # Note that double slashes can cause group-creation failures so we need
        # to carefully strip them out.
        bucket = (default_s3_path).rstrip("/")
        stamp = datetime.datetime.today().strftime("%Y%m%d-%H%M%S")
        creation_uri = f"tiledb://{namespace}/{bucket}/{dirname}/{basename}_{stamp}"
        readback_uri = f"tiledb://{namespace}/{basename}_{stamp}"
        return (creation_uri, readback_uri)

    uri = f"/tmp/tiledbsoma-cloud-test/{dirname}/{basename}"
    if pathlib.Path(uri).exists():
        shutil.rmtree(uri)
    pathlib.Path(pathlib.Path(uri).parent).mkdir(parents=True, exist_ok=True)
    # Please leave this comment in place.
    print()
    print("USING LOCAL URI", uri)
    print()
    return (uri, uri)


def util_tear_down_uri(uri):
    # This assumes tiledb.cloud.login has already been called at util_make_uri.
    if uri.startswith("tiledb://"):
        tiledb.cloud.groups.delete(uri=uri, recursive=True)
    # Delete local URIs only on _next_ run, so devs can inspect


def util_pbmc3k_unprocessed_versions():
    # New shape as in https://github.com/single-cell-data/TileDB-SOMA/issues/2407
    # which was released with tiledbsoma 1.15.0.
    return [
        ["tiledb://unittest/pbmc3k_unprocessed_1_7_3", {"shape": "old"}],
        ["tiledb://unittest/pbmc3k_unprocessed_1_12_3", {"shape": "old"}],
        ["tiledb://unittest/pbmc3k_unprocessed_1_14_5", {"shape": "old"}],
        ["tiledb://unittest/pbmc3k_unprocessed_1_15_0", {"shape": "new"}],
        ["tiledb://unittest/pbmc3k_unprocessed_1_15_7", {"shape": "new"}],
    ]
