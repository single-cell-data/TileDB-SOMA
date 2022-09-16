#!/bin/bash
set -euo pipefail
# Context: .github/workflows/*.yaml
apt-get update -y && apt-get install -y libgeos-dev
