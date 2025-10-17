# Carrara Tess

TileDB v3 aka Carrara introduces a new URL schema which supports relative paths and group membership lifecycle changes.

This directory contains unit tests for Carrara-specific behavior.

NB: this is EXPERIMENTAL and subject to change.

## Running the tests

Prerequisites:

- You have a user account in a Carrara deployment
- You have set up your tiledb profile so that you can log into the account with a profile name
- You have created a teamspace for testing use

With that, set the environment variables and install package dependencies:

```bash
pip install tiledb_client
export CARRARA_TEST_PROFILE="..."           # Profile name
export CARRARA_TEST_WORKSPACE="..."         # Workspace name (must match workspace for the profile)
export CARRARA_TEST_TEAMSPACE="..."         # Teamspace to use for tests
```

And _manually_ run the tests from the top-level repo directory:

```bash
pytest apis/python/remote_tests/carrara/ --carrara
```

Tests will be run within the specified workspace/teamspace, i.e., objects will be created at `tiledb://workspace/teamspace/...`
