# Test string array

```
import tiledb, numpy as np
# array of a-z
d = np.array([chr(c) * (c-96) for c in range(97,123)])
s = tiledb.schema_like(d)
tiledb.Array.create("/tmp/b1", s)
with tiledb.open("/tmp/b1", "w") as A:
    A[:] = d
```
