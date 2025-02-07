# How to run these tests

```
export TILEDB_REST_TOKEN="..."   # Get the token for the Saas `unittest` user
unsetTILEDB_REST_PAYER_NAMESPACE # If you have that set
```

As of 2025-02-07, use Python 3.9 to run UDF tests; otherwise they will be skipped.

```
python -m pytest path/to/this/directory
```

# Test-data setup

This is what was done for initial setup of these tests, and what should be done for future releases.

```
export TILEDB_REST_TOKEN="..." # Get the token for the Saas `unittest` user
export TILEDB_REST_PAYER_NAMESPACE=unittest
```

Here are source data you can find in the sandbox account `unittest` space:

```
s3://tiledb-unittest/soma-prod-test-data/h5ad/pbmc3k_unprocessed.h5ad
s3://tiledb-unittest/soma-prod-test-data/h5ad/pbmc3k_processed.h5ad
```

Local copy:

```
aws s3 cp s3://tiledb-unittest/soma-prod-test-data/h5ad .
```

Then use `tiledbsoma.io.from_h5ad` with the following sources and data:

* Preferr a bare Docker image
* Repeat for all desired TileDB-SOMA versions:
  * `pip install tiledbsoma==1.15.7` (or whichever version)
  * Ingest to `s3://tiledb-unittest/soma-prod-test-data/1.15.7/pbmc3k_unprocessed_1.15.7`
  * Register this in the cloud UI
    * Note: as of 2025-02-07 the cloud UI disallows `.` in group names so register with name like `1_15_7`.
    * Tracked at [sc-63068](https://app.shortcut.com/tiledb-inc/story/63068/allow-in-registration-paths)
  * Do not ingest directly to `tiledb://unittest/s3://tiledb-unittest/soma-prod-test-data/1.15.7/pbmc3k_unprocessed_1.15.7` since this will use today's version of core server-side, and what we want to really test is data written entirely by the pip-installed versions of tiledbsoma and core.
