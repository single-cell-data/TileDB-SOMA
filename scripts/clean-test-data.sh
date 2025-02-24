#!/usr/bin/env bash
#
# A script to download and untar test data
#

set -euo pipefail

# Change directory to the `test` folder.
cd "$(dirname "$0")/../data"
rm -rf soco
rm -rf example-visium-v2
