#!/usr/bin/env bash
#
# A script to download and untar test data
#

# Change directory to the `test` data folder.
cd "$(dirname "$0")/../test"

# Clean and un-tar saco dataset.
rm -rf soco && tar zxf soco.tgz
