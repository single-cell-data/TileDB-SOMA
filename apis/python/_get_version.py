#!/usr/bin/env python3
"""
Generate RELEASE-VERSION file for scikit-build-core fallback.
This script reads the version from version.py and writes it to RELEASE-VERSION file.
Note: With setuptools_scm, version is dynamically derived from git tags,
but RELEASE-VERSION serves as a fallback for cases where git is not available.
"""

import pathlib
import sys

# Add the current directory to the path so we can import version
sys.path.insert(0, str(pathlib.Path(__file__).parent))

try:
    import version

    version_string = version.get_version()
except Exception as e:
    print(f"Warning: Could not get version from version.py: {e}", file=sys.stderr)
    version_string = "0.0.0.dev0"

# Write to RELEASE-VERSION file (matches version.py convention)
version_file = pathlib.Path(__file__).parent / "RELEASE-VERSION"
version_file.write_text(version_string + "\n")

print(version_string)
