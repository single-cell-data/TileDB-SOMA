---
name: Bug report
about: Create a report to help us improve
title: "[BUG]"
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Provide a code example and any sample input data (e.g. an H5AD) as an attachment to reproduce this behavior.

**Expected behavior**
A clear and concise description of what you expected to happen.

**Versions (please complete the following information):**
 - TileDB-SOMA version:
 - Language and language version (e.g. Python 3.8, R 4.2.2):
 - OS (e.g. MacOS, Ubuntu Linux):

How to get the above information in Python:
```
import tiledbsoma, tiledb, sys 
tiledbsoma.__version__ 
tiledb.__version__ 
sys.version_info
```

How to get the above information in R:
```
library(tiledbsoma)
library(tiledb)
utils::packageVersion("tiledbsoma")
utils::packageVersion("tiledb")
R.Version()$version.string
```

**Additional context**
Add any other context about the problem here.
