---
title: Cheat sheet
---

As mentioned above, most of Quarto editing is markdown editing. GitHub-flavored markdown will get
you most of the way there -- `h1`/`h2`/etc with `#`, `##`, etc., bullet lists, and the like. Here we
list out, for handy reference, some of our most common TileDB idioms all collected in one place.

## Images

```
![_Dense and sparse arrays_](images/dense-array-sparse-array.png)
```

![_Dense and sparse arrays_](images/dense-array-sparse-array.png)

## Callouts

```
:::{.callout-note}
Note: here is how to make a note.
:::
```

:::{.callout-note}
Note: here is how to make a note.
:::

```
:::{.callout-tip}
Tip: here is how to make a tip.
:::
```

:::{.callout-tip}
Tip: here is how to make a tip.
:::

```
:::{.callout-warning}
Warning: here is how to make a warning.
:::
```

:::{.callout-warning}
Warning: here is how to make a warning.
:::

```
:::{.callout-caution}
Caution: here is how to make a caution.
:::
```

:::{.callout-caution}
Caution: here is how to make a caution.
:::

```
:::{.callout-important}
Important: here is how to make an important.
:::
```

:::{.callout-important}
Important: here is how to make an important.
:::

## Text formatting

```
* **Bold:** This is bolded.
* *Italic:* This is italicized.
* `Code:` This is in-line code.
* [This is a hyperlink](https://tiledb.com)
```

* **Bold:** This is bolded.
* *Italic:* This is italicized.
* `Code:` This is in-line code.
* [This is a hyperlink](https://tiledb.com)

## Code blocks

<pre>
```C
#include <stdio.h>
#include <stdlib.h>
#include <tiledb/tiledb.h>

// Name of array.
const char* array_name = "example_array";
```
</pre>

```C
#include <stdio.h>
#include <stdlib.h>
#include <tiledb/tiledb.h>

// Name of array.
const char* array_name = "example_array";
```

## Tabsets

<pre>
::: {.panel-tabset}
### Python
```python
# This is a multi-line code snippet
import tiledb, tiledb.cloud
```

### R
```R
# This is a multi-line code snippet
library(tiledb)
library(tiledbcloud)
```
:::
</pre>

::: {.panel-tabset}
### Python
```python
# This is a multi-line code snippet
import tiledb, tiledb.cloud
```

### R
```R
# This is a multi-line code snippet
library(tiledb)
library(tiledbcloud)
```
:::

## More

Please see the [Quarto Guide](https://quarto.org/docs/authoring/markdown-basics.html).
