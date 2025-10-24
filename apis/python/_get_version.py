#!/usr/bin/env python3
"""Script to get version for scikit-build-core."""

import sys
from pathlib import Path

# Add current directory to path
sys.path.insert(0, str(Path(__file__).parent))

from version import get_version

if __name__ == "__main__":
    print(get_version())
