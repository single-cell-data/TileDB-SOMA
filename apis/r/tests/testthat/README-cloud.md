# TileDB Cloud Tests

These tests verify that the tiledbsoma R package works correctly with TileDB Cloud.

## Running the Tests

### Prerequisites

- You have a TileDB Cloud account
- You have set up an REST API token
- You have write permissions to a namespace with registered cloud credentials

### Environment Variables

Set the following environment variables:

```bash
export SOMA_TEST_CLOUD="true"                  # Required to enable cloud tests
export TILEDB_REST_TOKEN="..."                 # REST API token
export CLOUD_TEST_NAMESPACE="..."              # Namespace name (default: "TileDB-Inc")
export CLOUD_TEST_BUCKET="..."                 # S3 bucket with registered credentials
```

### Running Tests

Run the cloud tests from the repository root:

```bash
# Run only cloud tests
Rscript -e "testthat::test_local(path = 'apis/r', filter = 'cloud')"
```

## Helper Functions

The `helper-cloud.R` file provides utility functions for cloud tests:

- `skip_if_no_cloud()` - Skips tests if `SOMA_TEST_CLOUD` is not set to "true"
- `get_cloud_config()` - Returns configuration from environment variables
- `get_cloud_base_uri()` - Builds the base URI for test objects
- `cloud_path()` - Creates a unique path with automatic cleanup
- `cloud_unique_id()` - Generates a unique ID for test assets
- `get_test_seurat_object()` - Returns a simplified Seurat object for testing
