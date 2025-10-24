#!/bin/bash
# Test script for local builds with scikit-build-core and vcpkg
# This script tests the build process locally before CI

set -e

echo "=== Testing TileDB-SOMA Python build with scikit-build-core and vcpkg ==="

# Determine script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

echo "Script directory: $SCRIPT_DIR"
echo "Repository root: $REPO_ROOT"

# Step 1: Setup vcpkg
echo ""
echo "=== Step 1: Setting up vcpkg ==="
bash "$SCRIPT_DIR/setup_build_env.sh"

# Check if vcpkg was set up successfully
if [ -z "$CMAKE_TOOLCHAIN_FILE" ]; then
    if [ -f "$REPO_ROOT/vcpkg/scripts/buildsystems/vcpkg.cmake" ]; then
        export CMAKE_TOOLCHAIN_FILE="$REPO_ROOT/vcpkg/scripts/buildsystems/vcpkg.cmake"
        echo "Set CMAKE_TOOLCHAIN_FILE=$CMAKE_TOOLCHAIN_FILE"
    else
        echo "Warning: CMAKE_TOOLCHAIN_FILE not set and vcpkg not found"
    fi
fi

# Step 2: Generate VERSION file
echo ""
echo "=== Step 2: Generating VERSION file ==="
cd "$REPO_ROOT/apis/python"
python3 _get_version.py > VERSION
cat VERSION

# Step 3: Create a clean build environment
echo ""
echo "=== Step 3: Creating clean build environment ==="

# Remove any existing build artifacts
cd "$REPO_ROOT/apis/python"
rm -rf build dist *.egg-info _skbuild

# Step 4: Install build dependencies
echo ""
echo "=== Step 4: Installing build dependencies ==="
pip install --upgrade pip build scikit-build-core pybind11 cmake

# Step 5: Build in verbose mode
echo ""
echo "=== Step 5: Building package ==="
cd "$REPO_ROOT/apis/python"
python -m build --wheel -v

# Step 6: Install and test
echo ""
echo "=== Step 6: Installing and testing ==="
cd "$REPO_ROOT/apis/python"
pip install dist/*.whl --force-reinstall

# Step 7: Run smoke test
echo ""
echo "=== Step 7: Running smoke test ==="
python -c "import tiledbsoma; print('Import successful!'); print('Module file:', tiledbsoma.pytiledbsoma.__file__); tiledbsoma.show_package_versions()"

echo ""
echo "=== Build test completed successfully! ==="

