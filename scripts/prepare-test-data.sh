#!/usr/bin/env bash
#
# A script to download and extract test data. Skips data that is already present in the
# directory. Make sure to update `clean-test-data.sh` after modifying this script.
#
# See ../data/README.md for instructions on updating this script with new data.
#

set -euo pipefail

echo "Begin preparing data."

# Change directory to the `data` folder.
cd "$(dirname "$0")/../data"


# Extract saco dataset.
if [ -d ../data/soco ]; then
    echo "-- Skipping dataset 'data/soco'; directory 'data/soco' already exists."
else
    echo "-- Preparing dataset 'data/soco' ..."
    tar zxf ../test/soco.tgz
    echo "   ... finished preparing 'data/soco'."
fi


# Download and extract Visium v2 dataset.
name="example-visium-v2"
if [ -d $name ]; then
    echo "-- Skipping dataset 'data/$name'; directory 'data/$name' already exists."
else
    echo "-- Preparing dataset 'data/$name' ..."
    mkdir $name && cd $name
    wget https://github.com/single-cell-data/TileDB-SOMA-Test-Data/releases/download/dataset-2025-02-19/filtered_feature_bc_matrix.h5
    wget https://github.com/single-cell-data/TileDB-SOMA-Test-Data/releases/download/dataset-2025-02-19/raw_feature_bc_matrix.h5
    wget https://github.com/single-cell-data/TileDB-SOMA-Test-Data/releases/download/dataset-2025-02-19/spatial.tar.gz
    wget https://github.com/single-cell-data/TileDB-SOMA-Test-Data/releases/download/dataset-2025-02-19/filtered_visium2_loc.csv
    tar zxf spatial.tar.gz
    cd ..
    echo "   ... finished preparing dataset 'data/$name'."
fi


# Download and extract Visium v2 dataset.
name="example-visium-v1"
if [ -d $name ]; then
    echo "-- Skipping dataset 'data/$name'; directory 'data/$name' already exists."
else
    echo "-- Preparing dataset 'data/$name' ..."
    mkdir $name && cd $name
    wget https://github.com/single-cell-data/TileDB-SOMA-Test-Data/releases/download/dataset-2025-02-21/filtered_feature_bc_matrix.h5
    wget https://github.com/single-cell-data/TileDB-SOMA-Test-Data/releases/download/dataset-2025-02-21/raw_feature_bc_matrix.h5
    wget https://github.com/single-cell-data/TileDB-SOMA-Test-Data/releases/download/dataset-2025-02-21/filtered_visium1_loc.csv
    wget https://github.com/single-cell-data/TileDB-SOMA-Test-Data/releases/download/dataset-2025-02-21/spatial.tar.gz
    tar zxf spatial.tar.gz
    cd ..
    echo "   ... finished preparing dataset 'data/$name'."
fi

# Download and extra soma-experiment-versions
# Nominal use case in CI:
# * Download the .tgz
# * Untar it
# * Test cases use it
# For development:
# * If the directory is there, use it
# * Else if the .tgz is there, untar it to make the directory
# * Else get the .tgz from cloud storage
# * This enables the developer to try out any modified data files _before_ making
#   a new release tag on the TileDB-SOMA-Test-Data repo
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
