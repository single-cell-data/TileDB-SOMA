# TileDB-Cloud Tests

This directory contains unit tests against TileDB-CLoud.

NB: this is EXPERIMENTAL and subject to change.

## Running the tests

Prerequisites:

- You have a user account on a TileDB-Cloud deployment
- You have set up your tiledb profile so that you can log into the account with a profile name
- You have created a teamspace for testing use

With that, set the environment variables and install package dependencies:

```bash
pip install tiledb_client
export TILEDB_CLOUD_TEST_PROFILE="..."           # Profile name
```

And _manually_ run the tests from the top-level repo directory:

```bash
pytest apis/python/remote_tests/carrara/ --carrara
```

Tests will be run within the specified workspace/teamspace, i.e., objects will be created at `tiledb://workspace/teamspace/...`
