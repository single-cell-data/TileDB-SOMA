# Simple test data 1

For ManagedQuery basic functionality tests.

Creation of array:

- TileDB-Py 0.12.2
- Steps:

```
import tiledb, numpy as np
tiledb.from_numpy("simple/dim-uint64_attr-int64_10cells", np.arange(10))
```
