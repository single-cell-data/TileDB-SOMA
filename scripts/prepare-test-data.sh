#!/usr/bin/env bash
#
# A script to download and extract test data.Skips data that is already present in the
# directory.
#
# See ../data/README.md for instructions on updating this script with new data.
#

set -euo pipefail

echo "Begin preparing data."

# Change directory to the `data` folder.
cd "$(dirname "$0")/../data"


# Extract saco dataset.
if [ -d ../test/soco ]; then
    echo "-- Skipping dataset 'data/soco'; directory 'data/soco' already exists."
else
    echo "-- Preparing dataset 'data/soco' ..."
    tar zxf ../test/soco.tgz
    echo "   ... finished preparing 'test/soco.tgz'."
fi


# Download and extract Visium v2 dataset.
if [ -d example-visium-v2 ]; then
    echo "-- Skipping dataset 'data/example-visium-vs'; directory 'data/example-visium-v2' already exists."
else
    echo "-- Preparing dataset 'data/example-visium-v2' ..."
    mkdir example-visium-v2 && cd example-visium-v2
    wget https://github.com/single-cell-data/TileDB-SOMA-Test-Data/releases/download/dataset-2025-02-19/filtered_feature_bc_matrix.h5
    wget https://github.com/single-cell-data/TileDB-SOMA-Test-Data/releases/download/dataset-2025-02-19/raw_feature_bc_matrix.h5
    wget https://github.com/single-cell-data/TileDB-SOMA-Test-Data/releases/download/dataset-2025-02-19/spatial.tar.gz
    tar zxf spatial.tar.gz
    cd ..
    echo "   ... finished preparing dataset 'data/example-visium-v2'."
fi
