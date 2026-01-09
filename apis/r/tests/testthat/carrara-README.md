# Carrara Tests

TileDB v3 aka Carrara introduces a new URL schema which supports relative paths and group membership lifecycle changes.

This directory contains unit tests for Carrara-specific behavior in the R package.

## Test Files

- `test-carrara-01-create.R` - Tests for creating SOMA objects with Carrara URIs
- `test-carrara-02-delete.R` - Tests for deleting SOMA objects
- `test-carrara-03-uri-enforcement.R` - Tests for URI validation and enforcement
- `test-carrara-04-error-handling.R` - Tests for error handling with Carrara operations

## Running the Tests

### Prerequisites

- You have a user account in a Carrara deployment
- You have set up your tiledb profile so that you can log into the account with a profile name
- You have created a teamspace for testing use

### Environment Variables

Set the following environment variables:

```bash
export SOMA_TEST_CARRARA="true"                # Required to enable carrara tests
export CARRARA_TEST_PROFILE="..."              # Profile name (default: "carrara")
export CARRARA_TEST_WORKSPACE="..."            # Workspace name (default: "TileDB-Inc-Staging")
export CARRARA_TEST_TEAMSPACE="..."            # Teamspace to use for tests (default: "aaron-dev")
export CARRARA_TEST_FOLDER="..."               # Folder within teamspace (default: "remote_test")
export TILEDB_REST_SERVER_ADDRESS="..."        # REST server URL (default: "https://api.staging.tiledb.io")
```

### Running Tests

Run the carrara tests from the repository root:

```bash
# Run only carrara tests
Rscript -e "testthat::test_local(path = 'apis/r', filter = 'carrara')"
```

Tests will be run within the specified workspace/teamspace, i.e., objects will be created at `tiledb://workspace/teamspace/folder/...`

## Helper Functions

The `helper-carrara.R` file provides utility functions for carrara tests:

- `skip_if_no_carrara()` - Skips tests if `SOMA_TEST_CARRARA` is not set to "true"
- `get_carrara_config()` - Returns configuration from environment variables
- `with_carrara_env()` - Sets up carrara environment variables for a test scope
- `get_base_uri()` - Builds the base URI for test objects
- `carrara_array_path()` - Creates a unique array path with automatic cleanup
- `carrara_group_path()` - Creates a unique group path with automatic cleanup
