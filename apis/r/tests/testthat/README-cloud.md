# TileDB Cloud Tests

These tests verify that the tiledbsoma R package works correctly with TileDB Cloud (v2).

## Setup

### Prerequisites

- A TileDB Cloud account with a REST API token
- Write permissions to a namespace with registered cloud storage credentials (i.e., an S3 bucket)

### Create a TileDB Profile

Cloud tests authenticate using a [TileDB profile](../../dev/reports/tiledb-embedded-credential-and-profile-system.md) stored in `~/.tiledb/profiles.json`. Add a named profile with your Cloud server address and token:

```json
{
  "version": 1,
  "cloudStaging": {
    "rest.server_address": "https://api.dev.tiledb.io",
    "rest.token": "your-cloud-api-token"
  }
}
```

Replace `"cloudStaging"` with whatever name you prefer and update the server address and token for your target Cloud environment.

### Environment Variables

The following environment variables control test execution:

```bash
# Required
export SOMA_TEST_CLOUD="true"                       # Enable cloud tests
export CLOUD_TEST_BUCKET="s3://your-bucket/prefix"  # S3 bucket with registered credentials

# Optional
export CLOUD_TEST_PROFILE="cloudStaging"            # Profile name (default: "default")
export CLOUD_TEST_NAMESPACE="your-namespace"        # Namespace name (default: "TileDB-Inc")
```

**Note:** `TILEDB_REST_TOKEN` should *not* be set. The tests automatically unset it to avoid conflicts with the profile-based authentication.

## Running Tests

From the `apis/r` directory:

```bash
Rscript -e 'testthat::test_local(filter = "cloud")'
```

Or from the repository root:

```bash
Rscript -e 'testthat::test_local(path = "apis/r", filter = "cloud")'
```

## Helper Functions

The `helper-cloud.R` file provides utility functions for cloud tests:

- `skip_if_no_cloud()` - Skips tests unless `SOMA_TEST_CLOUD=true` and `CLOUD_TEST_BUCKET` is set
- `get_cloud_config()` - Returns configuration (profile, namespace, bucket) from environment variables
- `with_cloud_env()` - Sets `TILEDB_PROFILE_NAME` and unsets `TILEDB_REST_TOKEN` for the test scope
- `get_cloud_base_uri()` - Builds the base URI for test objects (`tiledb://<namespace>/<bucket>`)
- `cloud_path()` - Creates a unique cloud URI with automatic cleanup
- `get_test_seurat_object()` - Returns a simplified Seurat object for testing
