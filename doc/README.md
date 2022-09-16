# Basic on-laptop setup

```
pip install .
```

This is very important -- _for developers_, the nominal use-mode is `python setup.py develop` as documented in our [../apis/python/README.md](../apis/python/README.md). But this local build _will not find_ `tiledbsoma-py` from this local install. You must install `pip install apis/python so it can find Python source for document autogen, _and_ you must re-run `pip install apis/python after each and every source-file edit, even if you're just doing an edit-build-preview iteration loop in a sandbox checkout.

```
#!/bin/bash
set -euo pipefail
sphinx-build -E -T -b html -d foo/doctrees -D language=en doc/source doc/html
```

```
#!/bin/bash
set -euo pipefail
open doc/html/python-api.html
```
