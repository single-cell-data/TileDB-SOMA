#!/usr/bin/env bash
#
# A script to download and untar test data
#
# See [the README in the `data/` directory](../data/README.md) for instructions on
# updating this script with new data.
#

set -euo pipefail

# Change directory to the `data` folder.
cd "$(dirname "$0")/../data"

# Remove prepared test data.
rm -rf soco
rm -rf example-visium-v2
