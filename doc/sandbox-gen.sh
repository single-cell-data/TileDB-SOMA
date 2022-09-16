#!/bin/bash

# Please run from ..

set -euo pipefail

sphinx-build -E -T -b html -d doc/doctrees -D language=en doc/source doc/html
