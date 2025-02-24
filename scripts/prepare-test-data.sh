#!/usr/bin/env bash

#
# A script to download and extract test data.
#
# Skips data that is already present in the directory.
#

set -euo pipefail

echo "Begin preparing data."

# Change directory to the `test` folder.
cd "$(dirname "$0")/../test"

# Extract saco dataset.
if [ -d soco ]; then
    echo "-- Skipping dataset 'test/soco'; directory 'test/soco' already exists."
else
    echo "-- Preparing dataset 'test/soco' ..."
    tar zxf soco.tgz
    echo "   ... finished preparing 'test/soco.tgz'."
fi

# Change directory to the `data` folder.
if [ ! -d ../data ]; then
  echo "Could not find ../data from pwd$(pwd); exiting"
  exit 1
fi
cd ../data

# Download and extract
name="soma-experiment-versions"
echo "-- Preparing dataset 'data/$name' ..."
if [ -d $name ]; then
    echo "-- Skipping dataset 'data/$name'; directory 'data/$name' already exists."
else
    if [ ! -f $name.tgz ]; then
        wget https://github.com/single-cell-data/TileDB-SOMA-Test-Data/releases/download/dataset-2025-02-24/$name.tgz
    fi
    tar zxf $name.tgz
fi
echo "   ... finished preparing dataset 'data/$name'."
