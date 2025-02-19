#!/usr/bin/env bash
#
# A script to download and untar test data
#

# Change directory to the `test` folder.
cd "$(dirname "$0")/../test"

# Clean saco dataset.
rm -rf soco

# Change direcotry to the `data` folder.
cd ../data
rm -rf example-visium-v2
